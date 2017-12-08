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
#include "typhon.h"

#include <cassert>



namespace _TYPH_Internal {
namespace {

int
Get_MPI_Op(TYPH_Datatype type, TYPH_Op in, MPI_Op *out)
{
    // Select appropriate MPI_Op based on TYPH_Op and datatype
    switch (type) {
    case TYPH_DATATYPE_REAL:
    case TYPH_DATATYPE_INTEGER:
        switch (in) {
        case TYPH_OP_SUM:
            *out = MPI_SUM;
            break;

        case TYPH_OP_PROD:
            *out = MPI_PROD;
            break;

        case TYPH_OP_MAX:
            *out = MPI_MAX;
            break;

        case TYPH_OP_MIN:
            *out = MPI_MIN;
            break;

        default:
            *out = MPI_OP_NULL;
            return TYPH_FAIL;
        }
        break;

    case TYPH_DATATYPE_LOGICAL:
        switch (in) {
        case TYPH_OP_OR:
            *out = MPI_LOR;
            break;

        case TYPH_OP_XOR:
            *out = MPI_LXOR;
            break;

        case TYPH_OP_AND:
            *out = MPI_LAND;
            break;

        default:
            *out = MPI_OP_NULL;
            return TYPH_FAIL;
        }
        break;

    default:
        *out = MPI_OP_NULL;
        return TYPH_FAIL;
    }

    return TYPH_SUCCESS;
}



template <typename T>
int
Gather_0D(T in, T *out, MPI_Comm *comm)
{
    using namespace _TYPH_Internal;
    TYPH_Datatype const type = TYPH_Datatype_From_T(in);
    MPI_Datatype const mpi_type = MPI_Datatype_From_TYPH_Datatype(type);
    if (mpi_type == MPI_DATATYPE_NULL) return TYPH_FAIL;

    int mpi_err = MPI_Allgather(&in, 1, mpi_type, out, 1, mpi_type, *comm);
    return mpi_err == MPI_SUCCESS ? TYPH_SUCCESS : TYPH_FAIL;
}



template <typename T>
int
Gather_1D(T *in, int len, T *out, MPI_Comm *comm)
{
    using namespace _TYPH_Internal;
    TYPH_Datatype const type = TYPH_Datatype_From_T(*in);
    MPI_Datatype const mpi_type = MPI_Datatype_From_TYPH_Datatype(type);
    if (mpi_type == MPI_DATATYPE_NULL) return TYPH_FAIL;

    int mpi_err = MPI_Allgather(in, len, mpi_type, out, len, mpi_type, *comm);
    return mpi_err == MPI_SUCCESS ? TYPH_SUCCESS : TYPH_FAIL;
}



template <typename T>
int
Reduce_0D(T in, T *out, TYPH_Op op, MPI_Comm *comm)
{
    using namespace _TYPH_Internal;

    TYPH_Datatype const type = TYPH_Datatype_From_T(in);
    MPI_Datatype const mpi_type = MPI_Datatype_From_TYPH_Datatype(type);
    if (TYPH_ASSERT(
            mpi_type != MPI_DATATYPE_NULL,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Invalid MPI datatype") == TYPH_FAIL) {
        return TYPH_FAIL;
    }

    MPI_Op mpi_op;
    int typh_err = Get_MPI_Op(type, op, &mpi_op);
    if (TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Invalid MPI op") == TYPH_FAIL) {
        return TYPH_FAIL;
    }

    int mpi_err = MPI_Allreduce(&in, out, 1, mpi_type, mpi_op, *comm);
    if (TYPH_ASSERT(
            mpi_err == MPI_SUCCESS,
            ERR_MPI,
            TYPH_ERR_MPI,
            "MPI_Allreduce failed") == TYPH_FAIL) {
        return TYPH_FAIL;
    }

    return TYPH_SUCCESS;
}



template <typename T>
int
Reduce_1D(T *in, int len, T *out, TYPH_Op op, MPI_Comm *comm)
{
    using namespace _TYPH_Internal;
    TYPH_Datatype const type = TYPH_Datatype_From_T(*in);
    MPI_Datatype const mpi_type = MPI_Datatype_From_TYPH_Datatype(type);
    if (mpi_type == MPI_DATATYPE_NULL) return TYPH_FAIL;

    MPI_Op mpi_op;
    int irc = Get_MPI_Op(type, op, &mpi_op);

    int mpi_err = MPI_SUCCESS;
    if (irc == TYPH_SUCCESS) {
        mpi_err = MPI_Allreduce(in, out, len, mpi_type, mpi_op, *comm);
    }

    bool const success = (irc == TYPH_SUCCESS) && (mpi_err == MPI_SUCCESS);
    return success ? TYPH_SUCCESS : TYPH_FAIL;
}



template <typename T>
int
Reduce_2D(T *in, int const *dims, T *out, TYPH_Op op, MPI_Comm *comm)
{
    using namespace _TYPH_Internal;
    TYPH_Datatype const type = TYPH_Datatype_From_T(*in);
    MPI_Datatype const mpi_type = MPI_Datatype_From_TYPH_Datatype(type);
    if (mpi_type == MPI_DATATYPE_NULL) return TYPH_FAIL;

    MPI_Op mpi_op;
    int irc = Get_MPI_Op(type, op, &mpi_op);

    int mpi_err = MPI_SUCCESS;
    if (irc == TYPH_SUCCESS) {
        int const len = dims[0] * dims[1];
        mpi_err = MPI_Allreduce(in, out, len, mpi_type, mpi_op, *comm);
    }

    bool const success = (irc == TYPH_SUCCESS) && (mpi_err == MPI_SUCCESS);
    return success ? TYPH_SUCCESS : TYPH_FAIL;
}

} // namespace
} // namespace _TYPH_Internal

// XXX(timrlaw): Bit of a horrid way of doing it but saves on code-duplication
//               since templates can't be used in the public API
#define GEN_TYPH_GATHER(type) \
    int \
    TYPH_Gather(type *in, int const *dims, int rank, type *out) \
    { \
        using namespace _TYPH_Internal; \
        MPI_Comm *comm = &TYPH_CORE->Get_MPI_Runtime()->comm; \
        switch (rank) { \
        case 0: \
            return Gather_0D<type>(*in, out, comm); \
        case 1: \
            return Gather_1D<type>(in, *dims, out, comm); \
        default: \
            TYPH_ASSERT(false, ERR_USER, TYPH_ERR_INVALID_ARG, \
                     "Unsupported rank (" + std::to_string(rank) + ")"); \
        } \
        return TYPH_FAIL; \
    }

GEN_TYPH_GATHER(int)
GEN_TYPH_GATHER(double)
GEN_TYPH_GATHER(bool)

#undef GEN_TYPH_GATHER



#define GEN_TYPH_REDUCE(type) \
    int \
    TYPH_Reduce(type *in, int const *dims, int rank, type *out, TYPH_Op op) \
    { \
        using namespace _TYPH_Internal; \
        MPI_Comm *comm = &TYPH_CORE->Get_MPI_Runtime()->comm; \
        switch (rank) { \
        case 0: \
            return Reduce_0D<type>(*in, out, op, comm); \
        case 1: \
            return Reduce_1D<type>(in, *dims, out, op, comm); \
        case 2: \
            return Reduce_2D<type>(in, dims, out, op, comm); \
        default: \
            TYPH_ASSERT(false, ERR_USER, TYPH_ERR_INVALID_ARG, \
                     "Unsupported rank (" + std::to_string(rank) + ")"); \
        } \
        return TYPH_FAIL; \
    }

GEN_TYPH_REDUCE(int)
GEN_TYPH_REDUCE(double)
GEN_TYPH_REDUCE(bool)

#undef GEN_TYPH_REDUCE
