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
#include <vector>
#include <iostream>
#include <mpi.h>

#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "core.h"
#include "decomposition.h"



namespace _TYPH_Internal {
namespace {

/**
 * \ingroup internal
 *
 * Store created partitioning information.
 */
std::vector<Partition_Info *> partitions;

} // namespace

/**
 * \param [in]  partition_id    partitioning info id
 * \param [out] partition_info  pointer to partitioning info
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
Get_Partition_Info(
        int partition_id,
        Partition_Info **partition_info)
{
    if ((decltype(partitions.size())) partition_id >= partitions.size()) {
        return TYPH_FAIL;
    }

    if (partition_info != nullptr) {
        *partition_info = partitions[partition_id];
    }

    return TYPH_SUCCESS;
}

} // namespace _TYPH_Internal

/**
 * \param [out] partition_id    ID for resulting structure
 * \param [in]  el_shape
 * \param [in]  num_layers      number of *ghost* layers (excluding real layer)
 * \param [in]  num_el_total
 * \param [in]  num_nd_total
 * \param [in]  el_to_proc
 * \param [in]  nd_to_proc
 * \param [in]  el_loc_to_glob
 * \param [in]  nd_loc_to_glob
 * \param [in]  connectivity
 * \param [in]  name
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Set_Partition_Info(
        int *partition_id,
        TYPH_Shape el_shape,
        int num_layers,
        int const *num_el_total,
        int const *num_nd_total,
        int const *el_to_proc,
        int const *nd_to_proc,
        int const *el_loc_to_glob,
        int const *nd_loc_to_glob,
        int const *connectivity,
        char const *cname)
{
    using namespace _TYPH_Internal;

    auto constexpr IXe2p = &Index_2D<2>;
    auto constexpr IXn2p = &Index_2D<2>;

    std::string name;
    if (cname != nullptr) {
        name = std::string(cname);
    }

    int irc = TYPH_SUCCESS;
    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    Partition_Info *partition_info = new Partition_Info();
    partitions.push_back(partition_info);

    partition_info->partition_id = partitions.size() - 1;
    *partition_id = partition_info->partition_id;

    partition_info->name = name;

    int const num_el = num_el_total[num_layers];
    int const num_nd = num_nd_total[num_layers];

    // Set nodes per element, based on element shape
    switch (el_shape) {
    case TYPH_SHAPE_QUAD:
        partition_info->nodes_per_el = 4;
        break;

    default:
        TYPH_ASSERT(false, ERR_USER, TYPH_ERR_INVALID_ARG,
                "Unrecognised element shape");
    }

    partition_info->num_layers = num_layers;

    partition_info->num_el_total = new int[num_layers + 1];
    TYPH_ALLOC_CHECK(partition_info->num_el_total, "partition_info->num_el_total");
    std::copy(num_el_total, num_el_total + num_layers + 1,
            partition_info->num_el_total);

    partition_info->num_nd_total = new int[num_layers + 1];
    TYPH_ALLOC_CHECK(partition_info->num_nd_total, "partition_info->num_nd_total");
    std::copy(num_nd_total, num_nd_total + num_layers + 1,
            partition_info->num_nd_total);

    partition_info->el_loc_to_glob = new int[num_el];
    TYPH_ALLOC_CHECK(partition_info->el_loc_to_glob, "partition_info->el_loc_to_glob");
    std::copy(el_loc_to_glob, el_loc_to_glob + num_el,
            partition_info->el_loc_to_glob);

    partition_info->nd_loc_to_glob = new int[num_nd];
    TYPH_ALLOC_CHECK(partition_info->nd_loc_to_glob, "partition_info->nd_loc_to_glob");
    std::copy(nd_loc_to_glob, nd_loc_to_glob + num_nd,
            partition_info->nd_loc_to_glob);

    // el_to_proc not needed in serial runs
    if (mp->size > 1) {
        partition_info->el_to_proc = new int[2 * num_el];
        TYPH_ALLOC_CHECK(partition_info->el_to_proc, "partition_info->el_to_proc");
        std::copy(el_to_proc, el_to_proc + 2*num_el,
                partition_info->el_to_proc);
    }

    partition_info->nd_to_proc = new int[2 * num_nd];
    TYPH_ALLOC_CHECK(partition_info->nd_to_proc, "partition_info->nd_to_proc");
    std::copy(nd_to_proc, nd_to_proc + 2*num_nd, partition_info->nd_to_proc);

    if (connectivity != nullptr) {
        partition_info->connectivity =
            new int[partition_info->nodes_per_el * num_el];
        TYPH_ALLOC_CHECK(partition_info->connectivity, "partition_info->connectivity");
        std::copy(
                connectivity,
                connectivity + partition_info->nodes_per_el * num_el,
                partition_info->connectivity);
    }

    // Print warning if current processor owns all cells/nodes, including
    // ghosts! Could be due to slide lines on all proc boundaries so not
    // necessarily an error
    bool warn = true;
    for (int i = 0; i < num_el; i++) {
        if (partition_info->el_to_proc[IXe2p(0, i)] != mp->rank) {
            warn = false;
            break;
        }
    }

    if (warn) {
        TYPH_WARNING("Processor " + std::to_string(mp->rank) +
            " owns *all* cells in partition, including ghosts, or partition"
            " has no ghosts");
    }

    warn = true;
    for (int i = 0; i < num_nd; i++) {
        if (partition_info->nd_to_proc[IXn2p(0, i)] != mp->rank) {
            warn = false;
            break;
        }
    }

    if (warn) {
        TYPH_WARNING("Processor " + std::to_string(mp->rank) +
            " owns *all* nodes in partition, including ghosts, or partition"
            " has no ghosts");
    }

    return irc;
}
