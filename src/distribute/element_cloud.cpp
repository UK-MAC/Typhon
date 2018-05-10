/* @HEADER@
 * Crown Copyright 2018 AWE.
 *
 * This file is part of Typhon.
 *
 * Typhon is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 * 
 * Typhon is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along with
 * Typhon. If not, see http://www.gnu.org/licenses/.
 * @HEADER@ */
#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "distribute/element_cloud.h"

#include <algorithm>
#include <numeric>
#include <memory>
#include <fstream>

#include "distribute/displacements.h"
#include "distribute/layer_info.h"


//#include "utilities/utils/sort.h"



namespace _TYPH_Internal {

Element_Cloud::Element_Cloud(int nd_per_el) :
    len(5 + nd_per_el) // 5 is the number of fixed properties
{
}



int &
Element_Cloud::operator()(int i, int j)
{
    assert(0 <= i && i < len);
    assert(0 <= j);
    return data[Index_2D(i, j, len)];
}



int const &
Element_Cloud::operator()(int i, int j) const
{
    assert(0 <= i && i < len);
    assert(0 <= j);
    return data[Index_2D(i, j, len)];
}



void
Element_Cloud::Set_Dims(int glob_size, int nproc, int rank)
{
    // Compute cloud dimensions
    size = ((glob_size - 1) / nproc) + 1;
    my_offset = rank * size;
    my_size = std::max(std::min(glob_size - my_offset, size), 0);
}



void
Element_Cloud::Set_Data(int nproc, int rank, Layer_Info &layer_info,
        MPI_Comm *comm, int nel_glob, Displacements &sdispls,
        Displacements &rdispls)
{
    // Set element cloud data from layer_info.conn_data

    // Determine cloud dimensions
    Set_Dims(nel_glob, nproc, rank);

    // Determine how much data is to be sent to each PE in cloud
    sdispls.zeroCounts();
    for (int iel = 0; iel < layer_info.nel; iel++) {
        int const iel_glob = layer_info(Element_Cloud::CGEL, iel);
        int const ip = iel_glob / size;
        assert(ip >= 0 && ip < nproc);
        sdispls.counts[ip]++;
    }

    // Data is already ordered correctly as layer_info.conn_data is sorted in
    // global element number

    // Determine receive counts
    int const recv_count = exchangeCounts(sdispls, rdispls, comm);
    int *rdata = new int[len * recv_count];

    // Send data
    sdispls.multiplyCounts(len);
    rdispls.multiplyCounts(len);
    sdispls.update();
    rdispls.update();
    exchangeData(layer_info.conn_data, rdata, sdispls, rdispls, comm);

    // Allocate EC data to my size
    data = new int[len * my_size];

    // Set EC data from rdata
    std::fill(data, data + len * my_size, 0);
    int cnt = 0; // number of instance found where cell data has already been set
    #define IXec(i, j) (Index_2D(i, j, this->len))
    for (int i = 0; i < recv_count; i++) {
        int const j = rdata[IXec(Element_Cloud::CGEL, i)] - my_offset;

        if ((*this)(Element_Cloud::CGEL, j) != 0) cnt++;
        for (int k = 0; k < len; k++) {
            (*this)(k, j) = rdata[IXec(k, i)];
        }
    }
    #undef IXec
    delete[] rdata;

    // Sum cnt across PEs
    int total;
    int typh_err = TYPH_Reduce_i(&cnt, nullptr, 0, &total, TYPH_OP_SUM);
    if (typh_err != TYPH_SUCCESS) {
        assert(false && "unhandled error");
        return;
    }

    if (rank == 0 && total > 0) {
        // Error if any cell connectivity is duplicated, indicate invalid mesh
        // partitioning
        assert(false && "Elements have been found on more than one partition");
        return;
    }
}



/**
 * @fn      Element_Cloud::Redistribute
 * @brief   Ensure processors have all the correct elements in the cloud
 *
 * @param [in]      rank
 * @param [inout]   sdispls
 * @param [inout]   rdispls
 * @param [inout]   layer_info
 * @param [inout]   comm
 */
void
Element_Cloud::Redistribute(int rank, Displacements &sdispls,
        Displacements &rdispls, Layer_Info &layer_info, MPI_Comm *comm)
{
    auto constexpr IX2 = &Index_2D<2>;

    // Called if parent routine has to redistribute mesh
    // Get connectivity table from EC based on target PE in EC table (entry 3)

    // Deallocate current connectivity data
    layer_info.Deallocate_Connectivity();

    // EC_data(EC_COWNPE,:) are the target PEs. Determine how much data is to be
    // sent to each.
    sdispls.zeroCounts();
    for (int i = 0; i < my_size; i++) {
        if ((*this)(Element_Cloud::CGEL, i) >= 0) {
            int const ip = (*this)(Element_Cloud::COWNPE, i);
            sdispls.counts[ip]++;
        }
    }

    int const send_count = sdispls.sumCounts();

    // Fill and communicate sconn
    {
        // Allocate array for connectivity data to be sent
        std::unique_ptr<int[]> sconn(new int[len * send_count]);

        // Fill sconn
        #define IXec(i, j) (Index_2D(i, j, this->len))
        sdispls.update();
        for (int i = 0; i < my_size; i++) {
            if ((*this)(Element_Cloud::CGEL, i) >= 0) {
                int const ip = (*this)(Element_Cloud::COWNPE, i);

                for (int k = 0; k < len; k++) {
                    sconn[IXec(k, sdispls.displs[ip])] = (*this)(k, i);
                }

                sdispls.displs[ip]++;
            }
        }
        #undef IXec

        // Determine receive counts
        layer_info.nel = exchangeCounts(sdispls, rdispls, comm);

        // Allocate new connectivity array
        layer_info.Allocate_Connectivity(len, layer_info.nel);

        // Send in sconn receive in linfo%conndata
        sdispls.multiplyCounts(len);
        rdispls.multiplyCounts(len);
        sdispls.update();
        rdispls.update();
        exchangeData(sconn.get(), layer_info.conn_data, sdispls, rdispls, comm);
    }

    // Sort layer_info.conn_data by global element number (column 1)
    Sort_0(layer_info.conn_data, layer_info.conn_dims, 2);

    // Determine local element numbers
    for (int i = 0; i < layer_info.nel; i++) {
        layer_info(Element_Cloud::COWNLEL, i) = i;
        layer_info(Element_Cloud::COWNPE,  i) = rank;
    }

    {
        // Allocate array sdata to send local element numbers to EC
        std::unique_ptr<int[]> sdata(new int[2 * layer_info.nel]);

        // Array of global, local element numbers, already sorted by destination
        // PE in EC
        sdispls.zeroCounts();
        for (int i = 0; i < layer_info.nel; i++) {
            int const j = layer_info(Element_Cloud::CGEL, i);

            sdata[IX2(0, i)] = j;   // Global element number
            sdata[IX2(1, i)] = i;   // Local element number

            int const ip = j / size;
            sdispls.counts[ip]++;
        }

        // Determine rcounts
        int const recv_count = exchangeCounts(sdispls, rdispls, comm);

        // Allocate receive array rdata
        std::unique_ptr<int[]> rdata(new int[2 * recv_count]);

        // Send in sdata, receive in rdata
        sdispls.multiplyCounts(2);
        rdispls.multiplyCounts(2);
        sdispls.update();
        rdispls.update();
        exchangeData(sdata.get(), rdata.get(), sdispls, rdispls, comm);

        // Loop over rdata entries and set local element numbers
        for (int i = 0; i < recv_count; i++) {
            int const j = rdata[IX2(0, i)] - my_offset;
            (*this)(Element_Cloud::COWNLEL, j) = rdata[IX2(1, i)];
        }
    }
}



//void
//Element_Cloud::dump(std::string filename, Error &err) const
//{
    //std::ofstream of(filename.c_str());
    //if (!of.is_open()) {
        //FAIL_WITH_LINE(err, "ERROR: couldn't open file " + filename);
        //return;
    //}

    //of << size << "\n";
    //of << my_size << "\n";
    //of << my_offset << "\n";
    //of << "\n";

    //if (data != nullptr) {
        //of << "data\n";
        //for (int i = 0; i < my_size; i++) {
            //for (int j = 0; j < Element_Cloud::LEN; j++) {
                //of << data[IXec(j, i)] << "\n";
            //}
            //of << "---\n";
        //}
        //of << "\n";
    //}

    //of.close();
//}

} // namespace _TYPH_Internal
