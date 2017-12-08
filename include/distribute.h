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
#ifndef TYPHON_DISTRIBUTE_H
#define TYPHON_DISTRIBUTE_H



// -----------------------------------------------------------------------------
// Public Typhon API - mesh distribution
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

/**
 * \brief Distribute the specified mesh amongst the current processors, given a
 *        partitioning.
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
        int **nd_owner);

/** @} */



#endif // TYPHON_DISTRIBUTE_H
