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
#ifndef TYPHON_UTILITIES_H
#define TYPHON_UTILITIES_H

#include <string>
#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <cassert>
#include <functional>

#include "typhon.h"
#include "types.h"



namespace _TYPH_Internal {

/**
 * \addtogroup internal
 *
 * @{
 */

// -----------------------------------------------------------------------------
// Sorting functions
// -----------------------------------------------------------------------------
// Remove duplicates and sort. nn is the new length
template <typename T>
inline void
Sort_Unique(T *ia, int n, int &nn)
{
    std::sort(ia, ia + n);
    auto last = std::unique(ia, ia + n);
    std::fill(last, ia + n, T());
    nn = last - ia;
}

void
Sort_0(int *iarray, int const *dims, int rank);

void
Sort_1(int *iarray, int const *dims, int rank, int const *icol, int len);

void
Sort_1(int *iarray, int const *dims, int rank, int const *icol, int len,
        int *indx, int indx_len __attribute__((unused)));



// -----------------------------------------------------------------------------
// Conversion between C, MPI and Typhon datatypes
// -----------------------------------------------------------------------------
inline MPI_Datatype
MPI_Datatype_From_TYPH_Datatype(TYPH_Datatype datatype)
{
    switch (datatype) {
    case TYPH_DATATYPE_INTEGER: return MPI_INT;
    case TYPH_DATATYPE_REAL:    return MPI_DOUBLE;
    case TYPH_DATATYPE_LOGICAL: return MPI_INT;
    default:
        return MPI_DATATYPE_NULL;
    }
}



// -----------------------------------------------------------------------------
// 2D array indexing functions
// -----------------------------------------------------------------------------
template <int W>
inline int
#ifndef TYPHON_DEBUG
constexpr
#endif
Index_2D(int j, int i) { // flipped to handle the col-major indexing
#ifdef TYPHON_DEBUG
    assert(i >= 0);
    assert(j >= 0);
    assert(j < W);
#endif
    return i * W + j;
}

inline int
#ifndef TYPHON_DEBUG
constexpr
#endif
Index_2D(int j, int i, int W) { // flipped to handle the col-major indexing
#ifdef TYPHON_DEBUG
    assert(i >= 0);
    assert(j >= 0);
    assert(j < W);
#endif
    return i * W + j;
}



// -----------------------------------------------------------------------------
// Debug output
// -----------------------------------------------------------------------------
/** \brief Execute a function on each rank, but serialised with barriers. */
int
Each_Rank(std::function<int(int)> f);



// -----------------------------------------------------------------------------
// Error checking functions
// -----------------------------------------------------------------------------
void _Warning(std::string routine, int line, std::string desc);

// Prints a warning to stderr.
#define TYPH_WARNING(desc) \
    _TYPH_Internal::_Warning( \
            __PRETTY_FUNCTION__, \
            __LINE__, \
            desc)

int _Assert(bool condition, std::string routine, int err_type,
        TYPH_Error error_code, std::string desc, int line = -1);

// Wrapper for _Assert that passes function name and line info. This should be
// used instead of using _Assert directly for this reason.
#define TYPH_ASSERT(cond, err_type, err_code, desc) \
    _TYPH_Internal::_Assert( \
            cond, \
            __PRETTY_FUNCTION__, \
            err_type, \
            err_code, \
            desc, \
            __LINE__)

// As TYPH_ASSERT, except will return TYPH_FAIL on failure of the asserted
// condition (when in release mode)
#define TYPH_ASSERT_RET(cond, err_type, err_code, desc) \
    if (_TYPH_Internal::_Assert( \
            cond, \
            __PRETTY_FUNCTION__, \
            err_type, \
            err_code, \
            desc, \
            __LINE__) == TYPH_FAIL) return TYPH_FAIL;

void _Alloc_Check(void *ptr, std::string var);

// Checks whether memory was allocated successfully.
#define TYPH_ALLOC_CHECK(ptr, name) \
    _TYPH_Internal::_Alloc_Check( \
            ptr, \
            name)

/** @} */

} // namespace _TYPH_Internal



#endif // TYPHON_UTILITIES_H
