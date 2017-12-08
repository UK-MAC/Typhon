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
#ifndef TYPHON_EXCHANGE_H
#define TYPHON_EXCHANGE_H



// -----------------------------------------------------------------------------
// Public Typhon API - exchange related functions
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

/** \brief Execute the given communication phase. */
int
TYPH_Exchange(
        int phase_id);

/** \brief Start the given communication phase. */
int
TYPH_Start_Exchange(
        int phase_id);

/** \brief Finalise the given communication phase. */
int
TYPH_Finish_Exchange(
        int phase_id);

/** @} */



#endif // TYPHON_EXCHANGE_H
