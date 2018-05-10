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
#ifndef TYPHON_DECOMPOSITION_H
#define TYPHON_DECOMPOSITION_H

#include <string>

#include "typhon.h"



namespace _TYPH_Internal {

/**
 * \addtogroup internal
 *
 * @{
 */

/** Store information regarding mesh partitioning. */
struct Partition_Info
{
    int partition_id;               // partition id
    int num_layers;                 // number of ghost layers
    int nodes_per_el;               // depends on calling code

    int *num_el_total = nullptr;    // 0..num_layers
    int *num_nd_total = nullptr;    // 0..num_layers

    int *el_to_proc = nullptr;      // 2, num_el
    int *nd_to_proc = nullptr;      // 2, num_nd

    int *el_loc_to_glob = nullptr;  // num_el
    int *nd_loc_to_glob = nullptr;  // num_nd

    int *connectivity = nullptr;    // npe, num_el

    std::string name;
};

/** \brief Get information about a given partitioning. */
int
Get_Partition_Info(
        int partition_id,
        Partition_Info **partition_info = nullptr);

/** @} */

} // namespace _TYPH_Internal



#endif // TYPHON_DECOMPOSITION_H
