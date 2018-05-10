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
#include <cassert>

#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "core.h"



namespace _TYPH_Internal {

void
_Warning(std::string routine, int line, std::string desc)
{
    std::cerr << "[Typhon WARNING] in \"" << routine << "\" at line "
              << line << "\n";
    std::cerr << "\n";
    std::cerr << "\t" << desc << "\n";
    std::cerr << "\n";
}



int
_Assert(
        bool condition,
#ifndef TYPHON_DEBUG
        std::string routine __attribute__((unused)),
        int err_type __attribute__((unused)),
        TYPH_Error err_code __attribute__((unused)),
        std::string desc __attribute__((unused)),
        int line __attribute__((unused)))
#else
        std::string routine,
        int err_type,
        TYPH_Error err_code,
        std::string desc,
        int line)
#endif
{
#ifndef TYPHON_DEBUG
    return condition ? TYPH_SUCCESS : TYPH_FAIL;
#else
    if (condition) return TYPH_SUCCESS;

    // If here, then the check failed and there's an error
    std::cerr << "[Typhon ERROR] in \"" << routine << "\" at line "
              << line << "\n";
    std::cerr << "\n";
    std::cerr << "\t" << desc << "\n";
    std::cerr << "\n";
    std::cerr << "Error type: " << err_type << "\n";
    std::cerr << "Error code: " << err_code << "\n";
    std::cerr << "\n";
    std::cerr << "[Typhon INFO] Aborting...\n";

    int typh_err = TYPH_FAIL;
    if (TYPH_CORE != nullptr) {
        typh_err = TYPH_CORE->Comms_Abort(TYPH_FAIL);
    }

    std::exit(typh_err);
    return typh_err;
#endif
}



void
_Alloc_Check(void *ptr, std::string var)
{
    TYPH_ASSERT(
        ptr != nullptr,
        ERR_MEM,
        TYPH_ERR_MEMORY,
        "Failed to allocate memory for " + var);
}


// This sorts a 2D array. absolutely minging FIXME(timrlaw)
void
Sort_0(int *iarray, int const *dims, int rank)
{
    assert(rank == 2 && "not implemented other ranks");

    #define IX1(i) ((i)-1)
    #define IX(i, j) (Index_2D((i)-1, (j)-1, dims[0]))
    int *itemp = new int[dims[0]];

    int const n = dims[1];
    if (n < 2) return;

    int l = n / 2 + 1;
    int ir = n;
    int i, j;
    int z;

label10:
    if (l > 1) {
        l--;
        for (z = 1; z <= dims[0]; z++) {
            itemp[IX1(z)] = iarray[IX(z, l)];
        }

    } else {
        for (z = 1; z <= dims[0]; z++) {
            itemp[IX1(z)] = iarray[IX(z, ir)];
        }

        for (z = 1; z <= dims[0]; z++) {
            iarray[IX(z, ir)] = iarray[IX(z, 1)];
        }

        ir--;
        if (ir == 1) {
            for (z = 1; z <= dims[0]; z++) {
                iarray[IX(z, 1)] = itemp[IX1(z)];
            }
            return;
        }
    }

    i = l;
    j = 2 * l;

label20:
    if (j <= ir) {
        if (j < ir) {
            if (iarray[IX(1, j)] < iarray[IX(1, j+1)]) j++;
        }

        if (itemp[IX1(1)] < iarray[IX(1, j)]) {
            for (z = 1; z <= dims[0]; z++) {
                iarray[IX(z, i)] = iarray[IX(z, j)];
            }
            i = j;
            j *= 2;

        } else {
            j = ir + 1;
        }

        goto label20;
    }

    for (z = 1; z <= dims[0]; z++) {
        iarray[IX(z, i)] = itemp[IX1(z)];
    }
    goto label10;

    #undef IX
}



void
Sort_1(int *iarray, int const *dims, int rank, int const *icol, int len)
{
    assert(rank == 2 && "rank must be 2");

    #define IX1(i) ((i)-1)
    #define IX2(i, j) (Index_2D((i)-1, (j)-1, dims[0]))

    int *itemp = new int[dims[0]];
    int l, ir, i, j, k, n, m;
    bool z;

    n = dims[1];
    m = len;
    if (n < 2) return;
    l = n/2 + 1;
    ir = n;

    while (true) {
        if (l > 1) {
            l--;
            for (int p = 1; p <= dims[0]; p++) {
                itemp[IX1(p)] = iarray[IX2(p, l)];
            }

        } else {
            for (int p = 1; p <= dims[0]; p++) {
                itemp[IX1(p)] = iarray[IX2(p, ir)];
            }

            for (int p = 1; p <= dims[0]; p++) {
                iarray[IX2(p, ir)] = iarray[IX2(p, 1)];
            }

            ir--;
            if (ir == 1) {
                for (int p = 1; p <= dims[0]; p++) {
                    iarray[IX2(p, 1)] = itemp[IX1(p)];
                }
                return;
            }
        }

        i = l;
        j = 2 * l;

        while (true) {
            if (j > ir) break;
            if (j < ir) {
                z = false;
                for (k = 1; k <= m; k++) {
                    if (icol[IX1(k)] > 0) {
                        if (iarray[IX2(icol[IX1(k)], j  )] ==
                            iarray[IX2(icol[IX1(k)], j+1)]) continue;
                        if (iarray[IX2(icol[IX1(k)], j  )] >
                            iarray[IX2(icol[IX1(k)], j+1)]) break;

                    } else {
                        if (abs(iarray[IX2(-icol[IX1(k)], j  )]) ==
                            abs(iarray[IX2(-icol[IX1(k)], j+1)])) continue;
                        if (abs(iarray[IX2(-icol[IX1(k)], j  )]) >
                            abs(iarray[IX2(-icol[IX1(k)], j+1)])) break;
                    }

                    z = true;
                    break;
                }
                if (z) j++;
            }

            z = false;
            for (k = 1; k <= m; k++) {
                if (icol[IX1(k)] > 0) {
                    if (itemp[IX1(icol[IX1(k)])] ==
                        iarray[IX2(icol[IX1(k)], j)]) continue;
                    if (itemp[IX1(icol[IX1(k)])] >
                        iarray[IX2(icol[IX1(k)], j)]) break;

                } else {
                    if (abs(itemp[IX1(-icol[IX1(k)])]) ==
                        abs(iarray[IX2(-icol[IX1(k)], j)])) continue;
                    if (abs(itemp[IX1(-icol[IX1(k)])]) >
                        abs(iarray[IX2(-icol[IX1(k)], j)])) break;
                }

                z = true;
                break;
            }

            if (z) {
                for (int p = 1; p <= dims[0]; p++) {
                    iarray[IX2(p, i)] = iarray[IX2(p, j)];
                }
                i = j;
                j *= 2;

            } else {
                j = ir + 1;
            }

        } // inner while (true)

        for (int p = 1; p <= dims[0]; p++) {
            iarray[IX2(p, i)] = itemp[IX1(p)];
        }

    } // outer while (true)

    delete[] itemp;

    #undef IX2
    #undef IX1
}



void
Sort_1(int *iarray, int const *dims, int rank, int const *icol, int len,
        int *indx, int indx_len __attribute__((unused)))
{
    assert(rank == 2 && "rank must be 2");

    #define IX1(i) ((i)-1)
    #define IX2(i, j) (Index_2D((i) - 1, (j) - 1, dims[0]))

    int *itemp = new int[dims[0]];
    int l, ir, i, j, k, it, n, m;
    bool z;

    n = dims[1];
    m = len;
    if (n < 2) return;
    l = n/2 + 1;
    ir = n;

    while (true) {
        if (l > 1) {
            l--;
            for (int p = 1; p <= dims[0]; p++) {
                itemp[IX1(p)] = iarray[IX2(p, l)];
            }
            it = indx[IX1(l)];

        } else {
            for (int p = 1; p <= dims[0]; p++) {
                itemp[IX1(p)] = iarray[IX2(p, ir)];
            }

            it = indx[IX1(ir)];
            for (int p = 1; p <= dims[0]; p++) {
                iarray[IX2(p, ir)] = iarray[IX2(p, 1)];
            }
            indx[IX1(ir)] = indx[IX1(1)];

            ir--;
            if (ir == 1) {
                for (int p = 1; p <= dims[0]; p++) {
                    iarray[IX2(p, 1)] = itemp[IX1(p)];
                }
                indx[IX1(1)] = it;
                return;
            }
        }

        i = l;
        j = 2 * l;

        while (true) {
            if (j > ir) break;
            if (j < ir) {
                z = false;
                for (k = 1; k <= m; k++) {
                    if (icol[IX1(k)] > 0) {
                        if (iarray[IX2(icol[IX1(k)], j  )] ==
                            iarray[IX2(icol[IX1(k)], j+1)]) continue;
                        if (iarray[IX2(icol[IX1(k)], j  )] >
                            iarray[IX2(icol[IX1(k)], j+1)]) break;

                    } else {
                        if (abs(iarray[IX2(-icol[IX1(k)], j  )]) ==
                            abs(iarray[IX2(-icol[IX1(k)], j+1)])) continue;
                        if (abs(iarray[IX2(-icol[IX1(k)], j  )]) >
                            abs(iarray[IX2(-icol[IX1(k)], j+1)])) break;
                    }

                    z = true;
                    break;
                }
                if (z) j++;
            }

            z = false;
            for (k = 1; k <= m; k++) {
                if (icol[IX1(k)] > 0) {
                    if (itemp[IX1(icol[IX1(k)])] ==
                        iarray[IX2(icol[IX1(k)], j)]) continue;
                    if (itemp[IX1(icol[IX1(k)])] >
                        iarray[IX2(icol[IX1(k)], j)]) break;

                } else {
                    if (abs(itemp[IX1(-icol[IX1(k)])]) ==
                        abs(iarray[IX2(-icol[IX1(k)], j)])) continue;
                    if (abs(itemp[IX1(-icol[IX1(k)])]) >
                        abs(iarray[IX2(-icol[IX1(k)], j)])) break;
                }

                z = true;
                break;
            }

            if (z) {
                for (int p = 1; p <= dims[0]; p++) {
                    iarray[IX2(p, i)] = iarray[IX2(p, j)];
                }
                indx[IX1(i)] = indx[IX1(j)];
                i = j;
                j *= 2;

            } else {
                j = ir + 1;
            }

        } // inner while (true)

        for (int p = 1; p <= dims[0]; p++) {
            iarray[IX2(p, i)] = itemp[IX1(p)];
        }
        indx[IX1(i)] = it;

    } // outer while (true)

    #undef IX2
    #undef IX1
}



int
Each_Rank(std::function<int(int)> f)
{
    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    int typh_err = TYPH_SUCCESS;
    for (int ip = 0; ip < mp->size; ip++) {
        if (ip == mp->rank) {
            typh_err = f(mp->rank);
        }

        typh_err = TYPH_Barrier();
    }

    return typh_err;
}

} // namespace _TYPH_Internal

/**
 *  The reason for this existing is to give C/Fortran client code something
 *  to free allocated memory with, a C++ application might just as well use
 *  delete[] directly.
 */
int
TYPH_Free(int *ptr)
{
    delete[] ptr;
    return TYPH_SUCCESS;
}



char const *
TYPH_To_String_Datatype(TYPH_Datatype datatype)
{
    switch (datatype) {
    case TYPH_DATATYPE_NULL:    return "TYPH_DATATYPE_NULL";
    case TYPH_DATATYPE_REAL:    return "TYPH_DATATYPE_REAL";
    case TYPH_DATATYPE_INTEGER: return "TYPH_DATATYPE_INTEGER";
    case TYPH_DATATYPE_LOGICAL: return "TYPH_DATATYPE_LOGICAL";
    default:                    return "";
    }
}



char const *
TYPH_To_String_Auxiliary(TYPH_Auxiliary aux)
{
    switch (aux) {
    case TYPH_AUXILIARY_NONE: return "TYPH_AUXILIARY_NONE";
    default:                  return "";
    }
}



char const *
TYPH_To_String_Auxiliary(TYPH_Centring centring)
{
    switch (centring) {
    case TYPH_CENTRING_NODE: return "TYPH_CENTRING_NODE";
    case TYPH_CENTRING_CELL: return "TYPH_CENTRING_CELL";
    default:                 return "";
    }
}
