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
#ifndef TYPHON_QUANTITY_H
#define TYPHON_QUANTITY_H



// -----------------------------------------------------------------------------
// Public Typhon API - quantity related functions
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

/** \brief Add a quantity to a communication phase. */
int
TYPH_Add_Quant_To_Phase(
        int phase_id,
        int quant_id,
        int recv_quant_id = -1,
        int key_set_id = -1,
        int ghosts_min = -1,
        int ghosts_max = -1);

/** \brief Set the base address for the given quantity. */
int
TYPH_Set_Quant_Address(
        int quant_id,
        void *data,
        int const *dims,
        int rank);

/** @} */



#endif // TYPHON_QUANTITY_H
