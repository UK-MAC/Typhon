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
#include "distribute/layer_info.h"

#include <cassert>
#include <numeric>
#include <fstream>
#include <iostream>
#include <memory>

#include "distribute/displacements.h"
#include "distribute/ns.h"
#include "distribute/node_cloud.h"
#include "distribute/element_cloud.h"



namespace _TYPH_Internal {
namespace {

auto constexpr IXnc = &Index_2D<Node_Cloud::LEN>;
auto constexpr IX2  = &Index_2D<2>;

#define IXec(i, j) (Index_2D(i, j, ec.len))

void
Get_Cloud_El(
        Layer_Info &layer_info,
        int nproc,
        int rank __attribute__((unused)),
        MPI_Comm *comm,
        Element_Cloud const &ec,
        Node_Cloud const &nc,
        Layer_Info const &prev_layer_info,
        Displacements &sdispls,
        Displacements &rdispls,
        NS *ns)
{
    // Get from NC the elements surrounding prev_layer_info->nd_glob. Does not
    // return element on the same processor or elements that have been
    // previously returned to this process.
    layer_info.lay       = prev_layer_info.lay + 1;
    layer_info.nd_per_el = prev_layer_info.nd_per_el;
    layer_info.ndoff     = prev_layer_info.ndoff + prev_layer_info.nnd;
    layer_info.eloff     = prev_layer_info.eloff + prev_layer_info.nel;

    // Loop over prev_layer_info->nd_glob and count how may got to each PE in NC
    sdispls.zeroCounts();
    for (int i = 0; i < prev_layer_info.nnd; i++) {
        int const ind_glob = prev_layer_info.nd_glob[i];  // Global node index
        int const ip = ind_glob / nc.size;
        assert(ip >= 0 && ip < nproc);
        sdispls.counts[ip]++;
    }

    // Get receive counts
    int recv_count = exchangeCounts(sdispls, rdispls, comm);
    int *rlist = new int[recv_count];
    assert(rlist != nullptr);

    // Send in prev_layer_info->nd_glob recvieve in rlist
    sdispls.update();
    rdispls.update();
    exchangeData(prev_layer_info.nd_glob, rlist, sdispls, rdispls, comm);

    // Loop over receive nodes by processor and count potential number of
    // elements to be sent
    int send_count = 0;
    int ii = 0;
    for (int ip = 0; ip < nproc; ip++) {
        for (int i = ii; i < ii + rdispls.counts[ip]; i++) {
            int const n = rlist[i] - nc.my_offset;

            // Loop over element entries for node
            for (int iii = nc.indx[n]; iii < nc.indx[n+1]; iii++) {
                int const ipp = nc.data[IXnc(Node_Cloud::COWNPE, iii)];
                if (ip != ipp) {
                    send_count++; // Not sending back to the source process
                }
            }
        }

        ii += rdispls.counts[ip];
    }

    // Allocate send list
    int *slist = new int[send_count];
    std::fill(slist, slist + send_count, 0);

    // Loop over receive nodes by processor
    sdispls.zeroCounts();
    send_count = 0;
    ii = 0;
    for (int ip = 0; ip < nproc; ip++) {
        for (int i = ii; i < ii + rdispls.counts[ip]; i++) {
            int const n = rlist[i] - nc.my_offset;

            // Loop over element entries for node
            for (int iii = nc.indx[n]; iii < nc.indx[n+1]; iii++) {
                int const ipp = nc.data[IXnc(Node_Cloud::COWNPE, iii)];
                if (ip != ipp) {

                    // Not sending back to the source process
                    if (ns->New(ip, nc.data_ref[iii])) {

                        // Element with reference NC_data_ref(iii) has not been
                        // sent to PE ip before, so must send to PE ip.
                        // comms_ns_new(ip,NC_data_ref(iii)) will now record
                        // element with reference NC_data_ref(iii) as having
                        // been sent to PE.
                        sdispls.counts[ip]++;
                        slist[send_count] = nc.data[IXnc(Node_Cloud::CGEL, iii)];
                        send_count++;
                    }
                }
            }
        }

        ii += rdispls.counts[ip];
    }

    delete[] rlist;
    rlist = nullptr;

    // Determine receive counts
    recv_count = exchangeCounts(sdispls, rdispls, comm);
    rlist = new int[recv_count];

    // Send in slist, receive in rlist
    sdispls.update();
    rdispls.update();
    exchangeData(slist, rlist, sdispls, rdispls, comm);

    // Deallocate send list slist
    delete[] slist;
    slist = nullptr;

    // Sort by global element number
    //utils::kernel::sort0(rlist, &recv_count, 1);
    std::sort(rlist, rlist + recv_count);

    // Determine list of unique new elements (not in prev_layer_info->el_glob).
    // Makes use of fact the at both rlist and prev_layer_info->el_glob) are in
    // increasing order.
    int const nel0 = prev_layer_info.nel;
    int cnt = 0;
    int prev = -1;
    int ilast = 0;

    // Loop over elements received
    ii = 0;
    for (int i = 0; i < recv_count; i++) {
        int const iel_glob = rlist[i]; // Global element number received
        int jj = 0;
        if (iel_glob != prev) {
            // Different element, look for element in prev_layer_info->el_glob
            for (ii = ilast; ii < nel0; ii++) {
                if (iel_glob <= prev_layer_info.el_glob[ii]) break;
            }

            if (ii < nel0) {
                if (iel_glob == prev_layer_info.el_glob[ii]) {
                    // iel_glob has been found as prev_layer_info->el_glob[ii]
                    ii++;

                } else {
                    // iel_glob is not in prev_layer_info->el_glob, count as new
                    // local element
                    jj = iel_glob;
                    cnt++;
                }

            } else {
                // Gone off the end of prev_layer_info->el_glob, iel_glob is not
                // in it. Count as new local element.
                jj = iel_glob;
                cnt++;
            }

            ilast = ii;
            prev = iel_glob;
        }

        // jj==iel_glob if iel_glob is a new element, otherwise jj==0: will
        // allow us to identify new nodes in next part
        rlist[i] = jj;
    }

    // Set number of new elements
    layer_info.nel = cnt;

    // Create list of the new element. Determine how many are sent to each PE in
    // the EC
    slist = new int[layer_info.nel];
    cnt = 0;
    sdispls.zeroCounts();
    for (int i = 0; i < recv_count; i++) {
        int const j = rlist[i];
        if (j > 0) {
            slist[cnt] = j;
            int const ip = j / ec.size;
            assert(ip >= 0 && ip < nproc);
            sdispls.counts[ip]++;
            cnt++;
        }
    }

    // Deallocate rlist
    delete[] rlist;
    rlist = nullptr;

    // Determine receive counts
    recv_count = exchangeCounts(sdispls, rdispls, comm);
    rlist = new int[recv_count];

    // Send in slist, receive in rlist
    sdispls.update();
    rdispls.update();
    exchangeData(slist, rlist, sdispls, rdispls, comm);

    // Deallocates slist
    delete[] slist;
    slist = nullptr;

    // Allocate array to send connectivity information
    int *sconn = new int[ec.len * recv_count];

    // Set connectivity information to be sent
    for (int i = 0; i < recv_count; i++) {
        for (int j = 0; j < ec.len; j++) {
            sconn[IXec(j, i)] = ec.data[IXec(j, rlist[i] - ec.my_offset)];
        }
    }

    // Deallocate rlist
    delete[] rlist;
    rlist = nullptr;

    // Allocate for final connectivity data
    layer_info.Allocate_Connectivity(ec.len, layer_info.nel);

    // Send back connectivity information (scaled reverse of slist/rlist
    // transfer above)
    sdispls.multiplyCounts(ec.len);
    rdispls.multiplyCounts(ec.len);

    sdispls.update();
    rdispls.update();
    exchangeData(sconn, layer_info.conn_data, rdispls, sdispls, comm);

    // Deallocate sconn send array
    delete[] sconn;
    sconn = nullptr;

    // Sort layer_info->conn_data by global element number
    Sort_0(layer_info.conn_data, layer_info.conn_dims, 2);
}



void
Populate_Element_Global_Indices(Layer_Info &layer_info)
{
    // Get global list of elements layer_info.el_glob from
    // layer_info.conn_data which is assumed to be sorted by increasing global
    // element number
    layer_info.nel = layer_info.conn_dims[1];
    layer_info.el_glob = new int[layer_info.nel];

    for (int i = 0; i < layer_info.nel; i++) {
       layer_info.el_glob[i] = layer_info(Element_Cloud::CGEL, i);
    }
}



void
Populate_Node_To_Element_Mapping(
        Layer_Info const &layer_info,
        int *&nd_to_el,
        int *nd_to_el_dims)
{
    auto constexpr IX2 = &Index_2D<2>;

    // Create node to element connectivity array from layer_info.conn_data
    std::unique_ptr<int[]> nodes(new int[layer_info.nd_per_el]);

    // Determine node to element table length n2ecnt
    int nd_to_el_size = 0;
    for (int i = 0; i < layer_info.nel; i++) {

        // Get list of nodes in element connectivity
        for (int j = 0; j < layer_info.nd_per_el; j++) {
            nodes[j] = layer_info(Element_Cloud::CCONS + j, i);
        }

        // There are cnt unique nodes in element connectivity
        int cnt;
        Sort_Unique(&nodes[0], layer_info.nd_per_el, cnt);

        // Add to total
        nd_to_el_size += cnt;
    }

    // Allocate node to element connectivity table
    nd_to_el_dims[0] = 2;
    nd_to_el_dims[1] = nd_to_el_size;

    nd_to_el = new int[nd_to_el_dims[0] * nd_to_el_dims[1]];

    // Loop over elements in connectivity table
    nd_to_el_size = 0;
    for (int i = 0; i < layer_info.nel; i++) {

        // Get list of nodes in element connectivity
        for (int j = 0; j < layer_info.nd_per_el; j++) {
            nodes[j] = layer_info(Element_Cloud::CCONS + j, i);
        }

        // Get the cnt unique nodes in element connectivity
        int cnt;
        Sort_Unique(&nodes[0], layer_info.nd_per_el, cnt);

        // Create entry for each of the cnt unique nodes
        for (int j = 0; j < cnt; j++) {
            int const n = nodes[j];

            nd_to_el[IX2(0, nd_to_el_size)] = n;
            nd_to_el[IX2(1, nd_to_el_size)] = i;
            nd_to_el_size++;
        }
    }

    // Sort by global node number
    Sort_0(nd_to_el, nd_to_el_dims, 2);
}



void
Get_New_Nodes(
        int *nd_to_el,
        int const *nd_to_el_dims,
        Layer_Info &layer_info,
        Layer_Info const *prev_layer_info)
{
    assert(nd_to_el != nullptr && "nd_to_el not provided");
    assert(nd_to_el_dims != nullptr && "nd_to_el_dims not provided");
    assert((nd_to_el_dims[0] == 2 || nd_to_el_dims[0] == Node_Cloud::LEN) &&
            "incorrectly sized nd_to_el");

    auto IXn2e = IX2;
    if (nd_to_el_dims[0] == Node_Cloud::LEN) {
        IXn2e = IXnc;
    }

    // Determine nodes in nd_to_el that are not in prev_layer_info->nd_glob and
    // create new local node numbers for these. Make use of the fact that both
    // nd_to_el and prev_layer_info.nd_glob are sorted. Modify the connectivity
    // layer_info.conn_data to be in terms of the new local numbering.

    // Can handle the case where prev_layer_info is not set
    int nnd0 = 0;
    if (prev_layer_info != nullptr) {
        nnd0 = prev_layer_info->nnd;
    }

    int prev = -1;
    int cnt = 0;
    int ilast = 0;

    int nn, iloc, j;
    for (int i = 0; i < nd_to_el_dims[1]; i++) {
        int const n = nd_to_el[IXn2e(0, i)];  // The global node number

        if (n != prev) { // This is a new node
            nn = n;

            // Look for node in prev_layer_info->nd_glob
            int ii = 0;
            for (ii = ilast; ii < nnd0; ii++) {
                if (n <= prev_layer_info->nd_glob[ii]) break;
            }

            if (ii < nnd0) {
                if (n == prev_layer_info->nd_glob[ii]) {
                    // Node found in prev_layer_info->nd_glob, set iloc to the
                    // local number of this node
                    iloc = prev_layer_info->ndoff + ii;
                    ii++;
                    nn = -n;

                } else {
                    // This is is a new node, create new local number
                    iloc = layer_info.ndoff + cnt;
                    cnt++;
                }

            } else {
                // Gone off the end of prev_layer_info->nd_glob, so this is is a
                // new node: create new local number
                iloc = layer_info.ndoff + cnt;
                cnt++;
            }

            ilast = ii; // Start point in prev_layer_info->nd_glob next time
            prev = n;
        }

        // This is the local element number
        j = nd_to_el[IXn2e(1, i)];

        // Replace any occurances of node n in connectivity with -iloc the local
        // node number
        if (j > -1) {
            for (int k = Element_Cloud::CCONS; k < layer_info.conn_dims[0]; k++) {
                if (layer_info(k, j) == n) {
                    layer_info(k, j) = -iloc;
                }
            }
        }

        // If this is not a new node then nn=-n, otherwise nn=n
        // This will enable us to easily identify the new nodes later
        nd_to_el[IXn2e(0, i)] = nn;
    }

    layer_info.nnd = cnt;
    layer_info.nd_glob = new int[layer_info.nnd];
    cnt = 0;
    prev = -1;
    for (int i = 0; i < nd_to_el_dims[1]; i++) {
        int const n = nd_to_el[IXn2e(0, i)];
        if ((n > -1) && (n != prev)) {
            layer_info.nd_glob[cnt++] = n;
            prev = n;
        }
    }

    cnt = 0;
    for (int i = 0; i < layer_info.conn_dims[1]; i++) {
        for (int k = Element_Cloud::CCONS; k < layer_info.conn_dims[0]; k++) {
            if (layer_info(k, j) > 0) cnt++;
        }
    }

    if (cnt > 0) {
        assert(false && "Local nodes not found in connectivity for layer");
        return;
    }

    // Negate connectivity
    for (int i = 0; i < layer_info.conn_dims[1]; i++) {
        for (int k = Element_Cloud::CCONS; k < layer_info.conn_dims[0]; k++) {
            layer_info(k, i) *= -1;
        }
    }
}



/**
 * @fn      Get_Base_Layer_Info
 * @brief   Generate base layer information from mesh connectivity and
 *          partitioning.
 *
 * @param [in]    nel_local  # elements in conn_data provided by this rank
 * @param [in]    nel_glob   global mesh element count
 * @param [in]    nnd_glob   global mesh node count
 * @param [in]    nd_per_el  nodes per element
 * @param [in]    nproc      MPI comm size
 * @param [in]    rank       MPI rank
 * @param [in]    comm       MPI communicator
 * @param [in]    conn_data  input connectivity data
 * @param [in]    partition  mesh partitioning
 * @param [inout] nc         node cloud
 * @param [inout] ec         element cloud
 * @param [inout] ns
 * @param [out]   layer_info layer info that we are generating
 */
void
Get_Base_Layer_Info(
        int nel_local,
        int nel_glob,
        int nnd_glob,
        int nd_per_el,
        int nproc,
        int rank,
        MPI_Comm *comm,
        int const *conn_data,
        int const *partition,
        Node_Cloud &nc,
        Element_Cloud &ec,
        NS *ns,
        Layer_Info &layer_info)
{
    Displacements sdispls(nproc);
    Displacements rdispls(nproc);

    int const nel = nel_local;

    // Initialise basic layer metadata
    layer_info.lay = 0;
    layer_info.nel = nel;
    layer_info.nd_per_el = nd_per_el;
    layer_info.Allocate_Connectivity((3 + nd_per_el) + 2, nel);

    // Copy connectivity data to layer info, leaving space for owner PEs and
    // local element numbering
    #define IXconn(i, j) (Index_2D((i), (j), (3 + nd_per_el)))
    int constexpr CONS = 3;
    for (int i = 0; i < nel; i++) {
        for (int j = 0; j < CONS; j++) {
            layer_info(j, i) = conn_data[IXconn(j, i)];
        }

        for (int j = CONS; j < 3 + nd_per_el; j++) {
            layer_info(j+2, i) = conn_data[IXconn(j, i)];
        }
    }
    #undef IXconn

    // Set owner PEs based on mesh partitioning, leave local element numbering
    // blank for now
    for (int i = 0; i < nel; i++) {
        layer_info(Element_Cloud::COWNPE, i) = partition[i];
        layer_info(Element_Cloud::COWNLEL, i) = -1;
    }

    // Sort connectivity by global element number
    Sort_0(layer_info.conn_data, layer_info.conn_dims, 2);

    // Send connectivity data to the element cloud
    ec.Set_Data(nproc, rank, layer_info, comm, nel_glob, sdispls, rdispls);

    // When redistributing mesh, replace layer_info.conn_data with connectivity
    // data for this PE. Create local element numbering and send back to EC
    ec.Redistribute(rank, sdispls, rdispls, layer_info, comm);

    // Get list of global elements
    Populate_Element_Global_Indices(layer_info);

    // Compute expanded node to element connectivity (nd_to_el) and send to
    // node cloud
    int *nd_to_el = nullptr;
    int nd_to_el_dims[2];
    nc.Set_Data(nproc, rank, comm, layer_info, nnd_glob, nd_to_el,
            nd_to_el_dims, ns);

    // Get list of global nodes and convert connectivity to local node numbering
    Get_New_Nodes(nd_to_el, nd_to_el_dims, layer_info, nullptr);

    // Deallocate node to element connectivity table
    delete[] nd_to_el;
}



/**
 * @fn      Get_Next_Layer_Info
 * @brief   Generate upper layer information from previous layer.
 *
 * @param [in]    nproc      MPI comm size
 * @param [in]    rank       MPI rank
 * @param [in]    comm       MPI communicator
 * @param [out]   layer_info previously generated layer info
 * @param [inout] nc         node cloud
 * @param [inout] ec         element cloud
 * @param [inout] ns         
 * @param [out]   layer_info layer info that we are generating
 */
void
Get_Next_Layer_Info(
        int nproc,
        int rank,
        MPI_Comm *comm,
        Layer_Info const &prev_layer_info,
        Node_Cloud &nc,
        Element_Cloud &ec,
        NS *ns,
        Layer_Info &layer_info)
{
    Displacements sdispls(nproc);
    Displacements rdispls(nproc);

    // Get connectivity for elements surrounding previous layer info
    Get_Cloud_El(layer_info, nproc, rank, comm, ec, nc, prev_layer_info,
            sdispls, rdispls, ns);

    // Get list of global elements layer_info.el_glob
    Populate_Element_Global_Indices(layer_info);

    // Get node to element connectivity (nd_to_el) for new elements
    int *nd_to_el = nullptr;
    int nd_to_el_dims[2];
    Populate_Node_To_Element_Mapping(layer_info, nd_to_el, nd_to_el_dims);

    // Get list of NEW global nodes layer_info.nd_glob and convert
    // layer_info.conn_data connectivity to local node numbering
    Get_New_Nodes(nd_to_el, nd_to_el_dims, layer_info, &prev_layer_info);

    delete[] nd_to_el;
}

} // namespace

Layer_Info::Layer_Info()
{
}



Layer_Info::~Layer_Info()
{
    if (el_glob != nullptr) delete[] el_glob;
    if (nd_glob != nullptr) delete[] nd_glob;

    Deallocate_Connectivity();
}



int &
Layer_Info::operator()(int i, int j)
{
    assert(0 <= i && i < conn_dims[0]);
    assert(0 <= j && j < conn_dims[1]);
    return conn_data[Index_2D(i, j, conn_dims[0])];
}



int const &
Layer_Info::operator()(int i, int j) const
{
    assert(0 <= i && i < conn_dims[0]);
    assert(0 <= j && j < conn_dims[1]);
    return conn_data[Index_2D(i, j, conn_dims[0])];
}



//void
//Layer_Info::Dump(std::string filename) const
//{
    //std::ofstream of(filename.c_str());
    //if (!of.is_open()) {
        //assert(false && "unhandled error");
        //return;
    //}

    //of << lay << "\n";
    //of << nel << "\n";
    //of << nnd << "\n";
    //of << eloff << "\n";
    //of << ndoff << "\n";
    //of << conn_dims[0] << "\n";
    //of << conn_dims[1] << "\n";
    //of << "\n";

    //if (el_glob != nullptr) {
        //of << "el_glob\n";
        //for (int i = 0; i < nel; i++) {
            //of << el_glob[i] << "\n";
        //}
        //of << "\n";
    //}

    //if (nd_glob != nullptr) {
        //of << "nd_glob\n";
        //for (int i = 0; i < nnd; i++) {
            //of << nd_glob[i] << "\n";
        //}
        //of << "\n";
    //}

    //if (conn_data != nullptr) {
        //of << "conn_data\n";
        //for (int i = 0; i < conn_dims[1]; i++) {
            //for (int j = 0; j < conn_dims[0]; j++) {
                //of << conn_data[Index_2D(j, i, conn_dims[0])] << "\n";
            //}
            //of << "---\n";
        //}
        //of << "\n";
    //}

    //of.close();
//}



void
Layer_Info::Allocate_Connectivity(int n1, int n2)
{
    Deallocate_Connectivity();

    conn_dims = new int[2];
    conn_dims[0] = n1;
    conn_dims[1] = n2;

    conn_data = new int[conn_dims[0] * conn_dims[1]];
}



void
Layer_Info::Deallocate_Connectivity()
{
    if (conn_data != nullptr) {
        delete[] conn_data;
        conn_data = nullptr;
    }

    if (conn_dims != nullptr) {
        delete[] conn_dims;
        conn_dims = nullptr;
    }
}



/**
 * @fn      Get_Layer_Info
 * @brief   Generate layer information from mesh connectivity and partitioning.
 *
 * @param [in]    nel_local      # elements in conn_data provided by this rank
 * @param [in]    nel_glob       global number of mesh elements
 * @param [in]    nnd_glob       global number of mesh nodes
 * @param [in]    nd_per_el      number of nodes per element
 * @param [in]    nproc          MPI comm size
 * @param [in]    rank           MPI rank
 * @param [in]    comm           MPI communicator
 * @param [in]    num_layers     how many layers we are generating
 * @param [in]    conn_data      input connectivity data
 * @param [in]    partition      mesh partitioning
 * @param [inout] nc             node cloud
 * @param [out]   layer_info_arr space for layer information
 */
void
Get_Layer_Info(
        int nel_local,
        int nel_glob,
        int nnd_glob,
        int nd_per_el,
        int nproc,
        int rank,
        MPI_Comm *comm,
        int num_layers,
        int const *conn_data,
        int const *partition,
        Node_Cloud &nc,
        Layer_Info *layer_info)
{
    Element_Cloud ec(nd_per_el);
    NS ns;

    // Build real base layer from connectivity data and mesh partitioning
    Get_Base_Layer_Info(
            nel_local,
            nel_glob,
            nnd_glob,
            nd_per_el,
            nproc,
            rank,
            comm,
            conn_data,
            partition,
            nc,
            ec,
            &ns,
            layer_info[0]);

    //TYPH_Barrier();
    //for (int ip = 0; ip < nproc; ip++) {
        //if (rank == ip) {
            //layer_info_arr[0].dump(
                    //"layer_info_arr[0]" + std::to_string(ip) + ".dat");
        //}
        //TYPH_Barrier();
    //}

    // Build each ghost layer, based off the layer below
    for (int l = 1; l <= num_layers; l++) {
        Get_Next_Layer_Info(
                nproc,
                rank,
                comm,
                layer_info[l-1],
                nc,
                ec,
                &ns,
                layer_info[l]);

        //TYPH_Barrier();
        //for (int ip = 0; ip < nproc; ip++) {
            //if (rank == ip) {
                //layer_info_arr[l].dump("layer_info_arr[" + std::to_string(l) +
                        //"]" + std::to_string(ip) + ".dat");
            //}
            //TYPH_Barrier();
        //}
    }

    //TYPH_Barrier();
    //for (int ip = 0; ip < nproc; ip++) {
        //if (rank == ip) {
            //nc.dump("node_cloud" + std::to_string(ip) + ".dat");
            //ec.dump("element_cloud" + std::to_string(ip) + ".dat");
        //}
        //TYPH_Barrier();
    //}

    // Cleanup data structures to track which elements have already been
    // returned by NC to this PE
    ns.Cleanup();
}

} // namespace _TYPH_Internal
