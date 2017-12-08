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
#ifndef TYPHON_KEYS_H
#define TYPHON_KEYS_H



namespace _TYPH_Internal {

/**
 * \addtogroup internal
 *
 * @{
 */

/** \brief Get a pointer to a created key set. */
int
Get_Key_Set(
        int key_set_id,
        Key_Set **key_set);

/** \brief Get a pointer to a key within a key set. */
int
Get_Key(
        Key_Set const *key_set,
        TYPH_Sendrecv sendrecv,
        int proc,
        int layer,
        int &key_id,
        Key **key);

/** @} */

} // namespace _TYPH_Internal

// -----------------------------------------------------------------------------
// Public Typhon API - key related functions
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

/** \brief Create a new key set. */
int
TYPH_Create_Key_Set(
        TYPH_Keytype key_type,
        int layer_min,
        int layer_max,
        int partition_id,
        int &key_set_id);

/** @} */



#endif // TYPHON_KEYS_H
