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
#include <string>
#include <vector>
#include <iostream>



namespace _TYPH_Internal {

Core *Core::singleton = nullptr;

void
Core::Init_Singleton()
{
    // Check if already initialised
    if (singleton != nullptr) return;

    singleton = new Core();
}



void
Core::Kill_Singleton()
{
    // Check if uninitialised
    if (singleton == nullptr) return;

    delete singleton;
}



int
Core::Initialise(
        MPI_Comm *comm)
{
    // Initialize the MPI communication
    return Comms_Init(comm);
}



int
Core::Kill(
        bool finalise)
{
    // Perform a barrier, so that all processes are okay and ready to exit
    // Typhon cleanly
    int irc = Comms_Barrier();

    // Call all registered kill routines for components - do in reverse order,
    // for cleanliness
    for (int i = kill_store.size() - 1; i >= 0; i--) {
        irc = kill_store[i]();
    }
    kill_store.clear();

    // Kill the MPI communication as needed, finalizing MPI if specified
    if (finalise) {
        irc = Comms_Kill(true);
    }

    mpi_runtime->initialised = false;
    return irc;
}



int
Core::Kill_Insert(
        Kill_Func kill_func)
{
    kill_store.push_back(kill_func);
    return TYPH_SUCCESS;
}



int
Core::Comms_Barrier()
{
    int mpi_err = MPI_Barrier(mpi_runtime->comm);
    return TYPH_ASSERT(
            mpi_err == MPI_SUCCESS,
            ERR_MPI,
            TYPH_ERR_MPI,
            "MPI_Barrier failure");
}



int
Core::Comms_Abort(
        int abort_code)
{
    mpi_runtime->finalised = true;

    // This is a terminal abort, so if supporting multiple communicators, make
    // sure they all get freed and then abort
    if (mpi_runtime->initialised) {
        MPI_Abort(mpi_runtime->comm, abort_code);
    }

    return 0;
}



Core::Core()
{
    mpi_runtime = new MPI_Runtime();
}



Core::~Core()
{
    delete mpi_runtime;
}



int
Core::Comms_Init(
        MPI_Comm *comm)
{
    // Initialize MPI communication, gets No. processes, task IDs, etc.
    int initialised;
    int mpi_err = MPI_Initialized(&initialised);
    TYPH_ASSERT_RET(
            mpi_err == MPI_SUCCESS,
            ERR_MPI,
            TYPH_ERR_MPI,
            "MPI_Initialized failure");

    if (!initialised) {
        mpi_err = MPI_Init(nullptr, nullptr);
        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Init failure");
    }

    MPI_Comm new_comm = MPI_COMM_WORLD;
    if (comm != nullptr) {
        new_comm = *comm;
    }

    int typh_err = Init_MPI_Runtime(&new_comm);
    return TYPH_ASSERT(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Init_MPI_Runtime failed");
}



int
Core::Comms_Kill(
        bool finalise)
{
    TYPH_ASSERT_RET(
            mpi_runtime->initialised,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "MPI_Runtime not initialised");

    if (finalise) {
        mpi_runtime->finalised = true;
        int mpi_err = MPI_Finalize();
        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Finalize failure");
    }

    return TYPH_SUCCESS;
}



int
Core::Init_MPI_Runtime(
        MPI_Comm *comm)
{
    if (*comm != MPI_COMM_NULL) {
        int mpi_err = MPI_Comm_dup(*comm, &mpi_runtime->comm);
        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Comm_dup failure");

        mpi_err = MPI_Comm_size(mpi_runtime->comm, &mpi_runtime->size);
        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Comm_size failure");

        mpi_err = MPI_Comm_rank(mpi_runtime->comm, &mpi_runtime->rank);
        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Comm_rank failure");

        mpi_runtime->min_rank    = 0;
        mpi_runtime->max_rank    = mpi_runtime->size - 1;
        mpi_runtime->master_rank = mpi_runtime->min_rank;
        mpi_runtime->is_master   = (mpi_runtime->rank == mpi_runtime->master_rank);

    } else {
        mpi_runtime->comm        = MPI_COMM_NULL;
        mpi_runtime->size        = 0;
        mpi_runtime->rank        = MPI_PROC_NULL;
        mpi_runtime->min_rank    = 0;
        mpi_runtime->max_rank    = 0;
        mpi_runtime->master_rank = 0;
        mpi_runtime->is_master   = false;
    }

    mpi_runtime->initialised = true;
    mpi_runtime->finalised   = false;

    return TYPH_SUCCESS;
}

} // namespace _TYPH_Internal

/**
 * \param [inout] comm  optionally use an existing MPI communicator
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Init(
        MPI_Comm *comm)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            !(TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Already initialised");

    Core::Init_Singleton();
    return Core::Get_Singleton()->Initialise(comm);
}



/**
 * \param [in] finalise whether to call MPI_Finalize
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Kill(
        bool finalise)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    int typh_err = TYPH_CORE->Kill(finalise);
    TYPH_ASSERT_RET(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Core::Kill failed");

    Core::Kill_Singleton();
    return TYPH_SUCCESS;
}



/**
 * \param [in] abort_code   exit code
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Abort(
        int abort_code)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    TYPH_CORE->Comms_Abort(abort_code);
    return TYPH_ASSERT(
            false,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Unreachable statement");
}



/**
 * \param [out] size    size of communicator
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Get_Size(
        int *size)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    TYPH_ASSERT_RET(
            size != nullptr,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Nullptr provided as argument");

    *size = TYPH_CORE->Get_MPI_Runtime()->size;
    return TYPH_SUCCESS;
}



/**
 * \param [out] rank    calling processors rank
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Get_Rank(
        int *rank)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    TYPH_ASSERT_RET(
            rank != nullptr,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Nullptr provided as argument");

    *rank = TYPH_CORE->Get_MPI_Runtime()->rank;
    return TYPH_SUCCESS;
}



/**
 * \returns if the calling processors is the master rank
 */
bool
TYPH_Is_Master()
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    return TYPH_CORE->Get_MPI_Runtime()->is_master;
}



/**
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Barrier()
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    return TYPH_CORE->Comms_Barrier();
}



/**
 * \returns time in seconds since some arbitrary point in the past
 */
double
TYPH_Get_Time()
{
    return MPI_Wtime();
}



/**
 * \param [out] comm    MPI communicator
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Set_Comm(
        MPI_Comm *comm)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            comm != nullptr,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Nullptr provided as argument");

    *comm = MPI_COMM_WORLD;
    return TYPH_SUCCESS;
}



/**
 * \param [out] comm    MPI self communicator
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Set_Comm_Self(
        MPI_Comm *comm)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            comm != nullptr,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Nullptr provided as argument");

    *comm = MPI_COMM_SELF;
    return TYPH_SUCCESS;
}
