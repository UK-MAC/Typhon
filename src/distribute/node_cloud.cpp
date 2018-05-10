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
#include "utilities.h"
#include "distribute/node_cloud.h"

#include <algorithm>
#include <memory>
#include <fstream>

#include "distribute/ns.h"
#include "distribute/displacements.h"
#include "distribute/element_cloud.h"
#include "distribute/layer_info.h"

//#include "dataAPI/params.h"
//#include "utilities/utils/sort.h"



namespace _TYPH_Internal {

auto constexpr IXnc = &Index_2D<Node_Cloud::LEN>;

Node_Cloud::~Node_Cloud()
{
    if (data != nullptr) {
        delete[] data;
        data = nullptr;
    }

    if (indx != nullptr) {
        delete[] indx;
        indx = nullptr;
    }

    if (data_ref != nullptr) {
        delete[] data_ref;
        data_ref = nullptr;
    }
}



void
Node_Cloud::Set_Dims(int glob_size, int nproc, int rank)
{
    // Compute cloud dimensions
    size = (glob_size - 1) / nproc + 1;
    my_offset = rank * (size);
    my_size = std::max(std::min(glob_size - my_offset, size), 0);
}



void
Node_Cloud::Set_Data(int nproc, int rank, MPI_Comm *comm,
        Layer_Info const &layer_info, int nnd_glob, int *&nd_to_el,
        int *nd_to_el_dims, NS *ns)
{
    Displacements sdispls(nproc);
    Displacements rdispls(nproc);

    // Construct extended node to element connectivity and send nodal data to
    // node cloud
    std::unique_ptr<int[]> nodes(new int[layer_info.nd_per_el]);

    // Determine node to element table length
    int nd_to_el_size = 0;
    for (int j = 0; j < layer_info.nel; j++) {

        // Get a list of nodes in element connectivity
        for (int k = 0; k < layer_info.nd_per_el; k++) {
            nodes[k] = layer_info(Element_Cloud::CCONS + k, j);
        }

        // There are cnt unique nodes in element connectivity
        int cnt;
        Sort_Unique(&nodes[0], layer_info.nd_per_el, cnt);

        // Add to total
        nd_to_el_size += cnt;
    }

    // Allocate node to element connectivity table
    nd_to_el_dims[0] = Node_Cloud::LEN;
    nd_to_el_dims[1] = nd_to_el_size;

    nd_to_el = new int[nd_to_el_dims[0] * nd_to_el_dims[1]];
    auto constexpr IXn2e = &Index_2D<Node_Cloud::LEN>;

    //n2ecnt=0_ink ! Total number of rows
    //n2ecnt2=0_ink ! Count number of positive node entries
    int num_rows = 0;           // Total number of rows
    int pve_nd_entries = 0;     // Count number of +ve node entries

    // Loop over elements in connectivity table
    for (int j = 0; j < layer_info.nel; j++) {

        // Get a list of nodes in element connectivity
        for (int k = 0; k < layer_info.nd_per_el; k++) {
            nodes[k] = layer_info(Element_Cloud::CCONS + k, j);
        }

        // Get the cnt unique nodes in element connectivity
        int cnt;
        Sort_Unique(&nodes[0], layer_info.nd_per_el, cnt);

        int const jg = layer_info(Element_Cloud::CGEL, j);

        // Create entry for each of the cnt unique nodes
        for (int k = 0; k < cnt; k++) {
            int const n = nodes[k];

            if (n >= 0) {
                pve_nd_entries++;

                nd_to_el[IXn2e(Node_Cloud::CGNOD, num_rows)] = n;
                nd_to_el[IXn2e(Node_Cloud::COWNLEL, num_rows)] = j;
            }

            nd_to_el[IXn2e(Node_Cloud::CGEL, num_rows)] = jg;
            nd_to_el[IXn2e(Node_Cloud::COWNPE, num_rows)] = rank;
            nd_to_el[IXn2e(Node_Cloud::COWNLNOD, num_rows)] = -1;

            num_rows++;
        }
    }

    // Sort by global node number
    Sort_0(nd_to_el, nd_to_el_dims, 2);

    // Set the local node numbers
    int prev = -1;
    int cnt = 0;
    for (int i = 0; i < num_rows; i++) {
        int const n = nd_to_el[IXn2e(Node_Cloud::CGNOD, i)];
        if (n != prev) {
            cnt++;
            prev = n;
        }

        nd_to_el[IXn2e(Node_Cloud::COWNLNOD, i)] = cnt - 1;
    }

    // Set node cloud dimensions
    Set_Dims(nnd_glob, nproc, rank);

    int *itmp = nullptr;
    if (pve_nd_entries < num_rows) {
        // There rows corresponding to negative node references. Transfer only
        // pve_nd_entries the positive rows of nd_to_el to new array itmp,
        // which is used instead of nd_to_el
        itmp = new int[Node_Cloud::LEN * pve_nd_entries];
        cnt = 0;
        for (int i = 0; i < num_rows; i++) {
            int const j = nd_to_el[IXn2e(Node_Cloud::COWNLEL, i)];
            if (j > 0) {
                for (int k = 0; k < Node_Cloud::LEN; k++) {
                    itmp[IXn2e(k, cnt)] = nd_to_el[IXn2e(k, i)];
                }

                cnt++;
            }
        }

        // Determine the number of entries to be sent to each PE in NC
        sdispls.zeroCounts();
        for (int i = 0; i < pve_nd_entries; i++) {
            int const n = itmp[IXn2e(Node_Cloud::CGNOD, i)];
            int const ip = n / size;
            assert(ip >= 0 && ip < nproc);
            sdispls.counts[ip]++;
        }

        // These are already stored in the correct order as sorted by global
        // node number

    } else {
        // Determine the number of entries to be sent to each PE in NC
        sdispls.zeroCounts();
        for (int i = 0; i < pve_nd_entries; i++) {
            int const n = nd_to_el[IXn2e(Node_Cloud::CGNOD, i)];
            int const ip = n / size;
            assert(ip >= 0 && ip < nproc);
            sdispls.counts[ip]++;
        }

        // These are already stored in the correct order as sorted by global
        // node number
    }

    // This is the amount of data to be stored in the node cloud on this PE
    my_data_size = exchangeCounts(sdispls, rdispls, comm);

    // Allocated data arrays
    data     = new int[Node_Cloud::LEN * my_data_size];
    data_ref = new int[my_data_size];
    std::fill(data_ref, data_ref + my_data_size, 0);

    // Send in nd_to_el, receive in data
    sdispls.multiplyCounts(Node_Cloud::LEN);
    rdispls.multiplyCounts(Node_Cloud::LEN);
    sdispls.update();
    rdispls.update();

    // Communicate either itmp or nd_to_el
    if (pve_nd_entries < num_rows) {
        exchangeData(itmp, data, sdispls, rdispls, comm);
        delete[] itmp;
    } else {
        exchangeData(nd_to_el, data, sdispls, rdispls, comm);
    }

    // Sort by global element number and set a unique reference number of
    // element received
    int dims[2];
    dims[0] = Node_Cloud::LEN;
    dims[1] = my_data_size;
    int col = Node_Cloud::CGEL + 1;
    Sort_1(data, dims, 2, &col, 1);

    prev = -1;
    cnt = 0;
    for (int i = 0; i < my_data_size; i++) {
        int const iel_glob = data[IXnc(Node_Cloud::CGEL, i)];
        if (iel_glob != prev) {
            prev = iel_glob;
            cnt++;
        }

        data_ref[i] = cnt - 1;
    }

    // cnt is the largest reference number, i.e. there are cnt unique elements
    // in the node cloud

    // Initialise NS* routines to store references from 1 to cnt. These routines
    // determine if a particular reference has been previously sent to a
    // particular PE number
    ns->Init(cnt);

    // Sort by global node number and global element number, applying the same
    // permutation to NC_data_ref
    int cols[2];
    cols[0] = Node_Cloud::CGNOD + 1;
    cols[1] = Node_Cloud::CGEL + 1;
    Sort_1(data, dims, 2, cols, 2, data_ref, my_data_size);

    // Create index of node data in node cloud on this PE
    indx = new int[my_size + 1];

    // Count number element entries for each node
    std::fill(indx, indx + my_size + 1, 0);
    for (int i = 0; i < my_data_size; i++) {
        int const n = data[IXnc(Node_Cloud::CGNOD, i)] - my_offset;
        indx[n]++;
    }

    // Construct index array (cumulative sum)
    int idx = 0;
    for (int i = 0; i < my_size; i++) {
        int const tmp = idx + indx[i];
        indx[i] = idx;
        idx = tmp;
    }
    indx[my_size] = idx;
}



//void
//Node_Cloud::dump(std::string filename, Error &err) const
//{
    //std::ofstream of(filename.c_str());
    //if (!of.is_open()) {
        //FAIL_WITH_LINE(err, "ERROR: couldn't open file " + filename);
        //return;
    //}

    //of << size << "\n";
    //of << my_size << "\n";
    //of << my_offset << "\n";
    //of << my_data_size << "\n";
    //of << "\n";

    //if (data != nullptr) {
        //of << "data\n";
        //for (int i = 0; i < my_data_size; i++) {
            //for (int j = 0; j < Node_Cloud::LEN; j++) {
                //of << data[IXnc(j, i)] << "\n";
            //}
            //of << "---\n";
        //}
        //of << "\n";
    //}

    //if (data_ref != nullptr) {
        //of << "data_ref\n";
        //for (int i = 0; i < my_data_size; i++) {
            //of << data_ref[i] << "\n";
        //}
        //of << "\n";
    //}

    //if (indx != nullptr) {
        //of << "indx\n";
        //for (int i = 0; i < my_size + 1; i++) {
            //of << indx[i] << "\n";
        //}
        //of << "\n";
    //}

    //of.close();
//}

} // namespace _TYPH_Internal
