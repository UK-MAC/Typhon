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
#ifndef TYPHON_COLLECT_H
#define TYPHON_COLLECT_H



// -----------------------------------------------------------------------------
// Public Typhon API - collect operations
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

/** Available reduction operations. */
enum TYPH_Op : int {
    TYPH_OP_SUM  = 1001,
    TYPH_OP_PROD = 1002,
    TYPH_OP_MAX  = 1003,
    TYPH_OP_MIN  = 1004,
    TYPH_OP_OR   = 1011,
    TYPH_OP_XOR  = 1012,
    TYPH_OP_AND  = 1013
};

/** Perform an all-gather operation on integer data. */
int
TYPH_Gather(
        int *in,
        int const *dims,
        int rank,
        int *out);

/** Perform an all-gather operation on real data. */
int
TYPH_Gather(
        double *in,
        int const *dims,
        int rank,
        double *out);

/** Perform an all-gather operation on boolean data. */
int
TYPH_Gather(
        bool *in,
        int const *dims,
        int rank,
        bool *out);

/** Perform an all-reduce operation on integer data. */
int
TYPH_Reduce(
        int *in,
        int const *dims,
        int rank,
        int *out,
        TYPH_Op op);

/** Perform an all-reduce operation on real data. */
int
TYPH_Reduce(
        double *in,
        int const *dims,
        int rank,
        double *out,
        TYPH_Op op);

/** Perform an all-reduce operation on boolean data. */
int
TYPH_Reduce(
        bool *in,
        int const *dims,
        int rank,
        bool *out,
        TYPH_Op op);

/** @} */



#endif // TYPHON_COLLECT_H
