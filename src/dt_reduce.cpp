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
#include <mpi.h>

#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "core.h"



namespace _TYPH_Internal {
namespace {

MPI_Datatype dt_reduce_type;
MPI_Op       dt_reduce_op;

// Procedure to perform the reduction operation: a min
void
Dt_Reduce_Op(
        void *in,
        void *in_out,
        int *len,
        MPI_Datatype *datatype __attribute__((unused)))
{
    TYPH_Dt *dt_in     = static_cast<TYPH_Dt *>(in);
    TYPH_Dt *dt_in_out = static_cast<TYPH_Dt *>(in_out);

    int const nsz = *len;
    for (int i = 0; i < nsz; i++) {
        if (dt_in_out[i].rdt > dt_in[i].rdt) {
            dt_in_out[i] = dt_in[i];
        }
    }
}

} // namespace
} // namespace _TYPH_Internal

// Procedure to register the type and operator with MPI
int
TYPH_Add_Reduce_Dt(void)
{
    using namespace _TYPH_Internal;

    int constexpr NTYPES = 4;

    int          block_len[NTYPES] = {0};
    MPI_Aint     disp[NTYPES] = {0};
    MPI_Datatype types[NTYPES] = {0};

    TYPH_Dt dummy;

    bool constexpr COMMUTATIVITY = true;

    int mpi_err = MPI_SUCCESS;
    mpi_err |= MPI_Get_address(&dummy.rdt, &disp[0]);
    mpi_err |= MPI_Get_address(&dummy.idt, &disp[1]);
    mpi_err |= MPI_Get_address(&dummy.sdt, &disp[2]);
    mpi_err |= MPI_Get_address(&dummy.mdt, &disp[3]);
    TYPH_ASSERT(mpi_err == MPI_SUCCESS, ERR_MPI, TYPH_ERR_MPI,
            "Failed to get MPI address displacements");

    for (int i = NTYPES - 1; i >= 0; i--) {
        disp[i] -= disp[0];
    }

    block_len[0] = 1;       // 1 double
    block_len[1] = 1;       // 1 integer
    block_len[2] = 8;       // 8 chars
    block_len[3] = 10;      // 10 chars

    types[0] = MPI_DOUBLE;
    types[1] = MPI_INT;
    types[2] = MPI_CHAR;
    types[3] = MPI_CHAR;

    // Create dt reduction datatype
    {
        mpi_err = MPI_Type_create_struct(
                NTYPES,
                block_len,
                disp,
                types,
                &_TYPH_Internal::dt_reduce_type);

        mpi_err |= MPI_Type_commit(
                &_TYPH_Internal::dt_reduce_type);

        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "Failed to create dt_reduce_type MPI derived type");
    }

    // Create dt reduction op
    {
        mpi_err = MPI_Op_create(
                &_TYPH_Internal::Dt_Reduce_Op,
                COMMUTATIVITY,
                &_TYPH_Internal::dt_reduce_op);

        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "Failed to create dt_reduce_op MPI op");
    }

    return TYPH_SUCCESS;
}



int
TYPH_Reduce_Dt(TYPH_Dt *val)
{
    using namespace _TYPH_Internal;

    TYPH_Dt rval;
    int mpi_err = MPI_Allreduce(val, &rval, 1, dt_reduce_type, dt_reduce_op,
            TYPH_CORE->Get_MPI_Runtime()->comm);

    TYPH_ASSERT_RET(
            mpi_err == MPI_SUCCESS,
            ERR_MPI,
            TYPH_ERR_MPI,
            "MPI_Allreduce failed");

    *val = rval;
    return TYPH_SUCCESS;
}
