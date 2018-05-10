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
#include "core.h"

#include "distribute/displacements.h"
#include "distribute/element_cloud.h"
#include "distribute/node_cloud.h"
#include "distribute/layer_info.h"



namespace _TYPH_Internal {
namespace {

/**
 * \fn      Get_Element_Owners
 * \brief   Extract element owner processors and local indices from layer info.
 *
 * \param [in]  layer_info  array of layer info structures for each ghost layer
 * \param [in]  num_layers  number of ghost layers
 * \param [out] el_owner    output array (2, nel)
 */
void
Get_Element_Owners(
        Layer_Info const *layer_info,
        int num_layers,
        int *el_owner)
{
    MPI_Runtime *mp = TYPH_CORE->Get_MPI_Runtime();

    auto constexpr IX2 = &Index_2D<2>;
    for (int l = 0; l <= num_layers; l++) {
        Layer_Info const &li = layer_info[l];

        int const eloff = li.eloff;
        for (int i = 0; i < li.nel; i++) {
            if (l == 0) {
                el_owner[IX2(0, eloff + i)] = mp->rank;
                el_owner[IX2(1, eloff + i)] = eloff + i;

            } else {
                el_owner[IX2(0, eloff + i)] = li(Element_Cloud::COWNPE, i);
                el_owner[IX2(1, eloff + i)] = li(Element_Cloud::COWNLEL, i);
            }
        }
    }
}



/**
 * \fn      Get_Node_Owners
 * \brief   Get node owner processors and local indices from node cloud.
 *
 * \param [in]  nc          node cloud
 * \param [in]  nd_loc_glob global indices for local nodes
 * \param [in]  total_nnd   total local nodes
 * \param [out] nd_owner    output array (2, nnd)
 */
void
Get_Node_Owners(
        Node_Cloud const &nc,
        int const *nd_loc_glob,
        int total_nnd,
        int *nd_owner)
{
    MPI_Runtime *mp = TYPH_CORE->Get_MPI_Runtime();

    // Construct sdata(1,:) a list of global nodes we need information for
    // Send each part to the appropriate PE in NC
    Displacements sdispls(mp->size);
    Displacements rdispls(mp->size);

    sdispls.zeroCounts();
    for (int i = 0; i < total_nnd; i++) {
        int const n = nd_loc_glob[i];
        int const ip = n / nc.size;
        assert(0 <= ip && ip < mp->size);
        sdispls.counts[ip]++;
    }
    sdispls.update();

    int *sdata = new int[2 * total_nnd];
    auto constexpr IX2  = &Index_2D<2>;

    for (int i = 0; i < total_nnd; i++) {
        int const n = nd_loc_glob[i];
        int const ip = n / nc.size;
        assert(0 <= ip && ip < mp->size);
        sdata[IX2(0, sdispls.displs[ip])] = n;
        sdata[IX2(1, sdispls.displs[ip])] = 0;
        sdispls.displs[ip]++;
    }

    int const recv_count = exchangeCounts(sdispls, rdispls, &mp->comm);
    int *rdata = new int[2 * recv_count];

    sdispls.multiplyCounts(2);
    rdispls.multiplyCounts(2);
    sdispls.update();
    rdispls.update();
    exchangeData(sdata, rdata, sdispls, rdispls, &mp->comm);

    // Loop over receive list rdata and set source PE and corresponding source
    // local node number for each entry
    auto constexpr IXnc = &Index_2D<Node_Cloud::LEN>;
    for (int i = 0; i < recv_count; i++) {
        int const n = rdata[IX2(0, i)] - nc.my_offset;

        // Last element surrounding node with have greatest global element
        // number as sorted by global node number and then global element number
        int const ii = nc.indx[n+1] - 1;
        rdata[IX2(0, i)] = nc.data[IXnc(Node_Cloud::COWNPE,   ii)];
        rdata[IX2(1, i)] = nc.data[IXnc(Node_Cloud::COWNLNOD, ii)];
    }

    // Send back data to requesting PE, into sdata
    exchangeData(rdata, sdata, rdispls, sdispls, &mp->comm);
    delete[] rdata;

    // Loop over sdata and set nd_owner values
    sdispls.divideCounts(2);
    sdispls.update();

    for (int i = 0; i < total_nnd; i++) {
        int const n = nd_loc_glob[i];
        int const ip = n / nc.size;
        assert(0 <= ip && ip < mp->size);

        nd_owner[IX2(0, i)] = sdata[IX2(0, sdispls.displs[ip])];
        nd_owner[IX2(1, i)] = sdata[IX2(1, sdispls.displs[ip])];
        sdispls.displs[ip]++;
    }

    delete[] sdata;
}

} // namespace
} // namespace _TYPH_Internal

/**
 * \param [in]  nel_glob    global mesh element count
 * \param [in]  nnd_glob    global mesh node count
 * \param [in]  nd_per_el   number of nodes per element
 * \param [in]  num_layers  number of ghost layers
 * \param [in]  conn_data   element connectivity
 * \param [in]  partition   owner processors for each global node
 * \param [out] layer_nel   number of elements per ghost layer
 * \param [out] layer_nnd   number of nodes per ghost layer
 * \param [out] total_nel   total number of elements across all layers
 * \param [out] total_nnd   total number of nodes across all layers
 * \param [out] el_loc_glob element global indices (total_nel)
 * \param [out] nd_loc_glob node global indices (total_nnd)
 * \param [out] el_region   element region indices (total_nel)
 * \param [out] el_material element material indices (total_nel)
 * \param [out] el_nd       element node global indices (nd_per_el, total_nel)
 * \param [out] el_owner    element owners and local indices (2, total_nel)
 * \param [out] nd_owner    node owners and local indices (2, total_nnd)
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Distribute_Mesh(
        int nel_local,
        int nel_glob,
        int nnd_glob,
        int nd_per_el,
        int num_layers,
        int const *conn_data,
        int const *partition,
        int *layer_nel,
        int *layer_nnd,
        int *total_nel,
        int *total_nnd,
        int **el_loc_glob,
        int **nd_loc_glob,
        int **el_region,
        int **el_material,
        int **el_nd,
        int **el_owner,
        int **nd_owner)
{
    using namespace _TYPH_Internal;

    TYPH_ASSERT_RET(
            nel_local > 0 && nel_glob > 0 && nnd_glob > 0,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Mesh must be non-zero size");

    TYPH_ASSERT_RET(
            nd_per_el >= 3,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Must be at least 3 nodes per element");

    TYPH_ASSERT_RET(
            num_layers >= 1,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Must be at least 1 ghost layer");

    MPI_Runtime *mp = TYPH_CORE->Get_MPI_Runtime();

    // Connectivity array conn_data entries:
    //
    //      - global element index
    //      - element region index
    //      - element material index
    //      - element node global index 1
    //      - element node global index 2
    //                  ...
    //      - element node global index n
    //

    // Get layer info
    Node_Cloud nc;
    Layer_Info *layer_info = new Layer_Info[num_layers + 1];
    Get_Layer_Info(
            nel_local,
            nel_glob,
            nnd_glob,
            nd_per_el,
            mp->size,
            mp->rank,
            &mp->comm,
            num_layers,
            conn_data,
            partition,
            nc,
            layer_info);

    // # local elements and nodes each layer
    for (int l = 0; l <= num_layers; l++) {
        layer_nel[l] = layer_info[l].nel;
        layer_nnd[l] = layer_info[l].nnd;
    }

    // Sum across all layers for totals
    *total_nel = std::accumulate(layer_nel, layer_nel + num_layers + 1, 0);
    *total_nnd = std::accumulate(layer_nnd, layer_nnd + num_layers + 1, 0);

    // Set element/node loc/glob mapping
    *el_loc_glob = new int[*total_nel];
    *nd_loc_glob = new int[*total_nnd];
    for (int l = 0; l <= num_layers; l++) {
        int const eloff = layer_info[l].eloff;
        for (int i = 0; i < layer_info[l].nel; i++) {
            (*el_loc_glob)[eloff + i] = layer_info[l].el_glob[i];
        }

        int const ndoff = layer_info[l].ndoff;
        for (int i = 0; i < layer_info[l].nnd; i++) {
            (*nd_loc_glob)[ndoff + i] = layer_info[l].nd_glob[i];
        }
    }

    // Set element region/material/nodes
    *el_region   = new int[*total_nel];
    *el_material = new int[*total_nel];
    *el_nd       = new int[nd_per_el * *total_nel];

    #define IXnpe(i, j) (Index_2D(i, j, nd_per_el))
    for (int l = 0; l <= num_layers; l++) {
        Layer_Info const &li = layer_info[l];

        int const eloff = layer_info[l].eloff;
        for (int i = 0; i < layer_info[l].nel; i++) {
            (*el_region)[eloff + i]   = li(Element_Cloud::CREG, i);
            (*el_material)[eloff + i] = li(Element_Cloud::CMAT, i);

            for (int k = 0; k < nd_per_el; k++) {
                int const kk = Element_Cloud::CCONS + k;
                (*el_nd)[IXnpe(k, eloff + i)] = li(kk, i);
            }
        }
    }
    #undef IXnpe

    // Set element owners/local indices
    *el_owner = new int[2 * *total_nel];
    Get_Element_Owners(layer_info, num_layers, *el_owner);

    // Set node owners/local indices
    *nd_owner = new int[2 * *total_nnd];
    Get_Node_Owners(nc, *nd_loc_glob, *total_nnd, *nd_owner);

    delete[] layer_info;
    return TYPH_SUCCESS;
}
