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
#ifndef TYPHON_SCHEDULE_H
#define TYPHON_SCHEDULE_H



namespace _TYPH_Internal {

/**
 * \addtogroup internal
 *
 * @{
 */

/** \brief Build the schedule for a given phase. */
int
Build_Schedule(
        int phase_id,
        int num_layers = -1);

/** \brief Delete the previously built schedule for a given phase. */
int
Delete_Schedule(
        int phase_id);

/** \brief Commit the previously built schedule for a given phase . */
int
Commit_Phase(
        Phase *phase);

/** @} */

} // namespace _TYPH_Internal



#endif // TYPHON_SCHEDULE_H
