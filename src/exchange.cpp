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
#include <mpi.h>

#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "core.h"
#include "register.h"
#include "schedule.h"



namespace _TYPH_Internal {
namespace {

/**
 * \param [in] phase_id     the phase to start the exchange for
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
Start_Exchange(int phase_id)
{
    int typh_err = TYPH_SUCCESS;
    MPI_Runtime const *mpi_runtime = TYPH_CORE->Get_MPI_Runtime();

    Phase *phase;
    typh_err = TYPH_REGISTRY->Get_Phase(phase_id, &phase);
    TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Phase " + std::to_string(phase_id) + " not found");

    // XXX this stuff doesn't appear to need to be done, although maybe better
    // to do this error checking up front?
    //Key_Set *key_set;
    //typh_err = Get_Key_Set(phase->key_set_id, &key_set);
    //TYPH_ASSERT(
            //typh_err == TYPH_SUCCESS,
            //ERR_INT,
            //TYPH_ERR_INTERNAL,
            //"Key set " + std::to_string(phase->key_set_id) + " not found");

    //Partition_Info *partition_info;
    //typh_err = Get_Partition_Info(key_set->partition_id, &partition_info);
    //TYPH_ASSERT(
            //typh_err == TYPH_SUCCESS,
            //ERR_INT,
            //TYPH_ERR_INTERNAL,
            //"Partition info " + std::to_string(key_set->partition_id) +
                //" not found");

    // Build the phase if required, should be done once for pure
    if (!phase->is_built) {
        if (phase->schedule != nullptr) {
            typh_err = Delete_Schedule(phase_id);
            TYPH_ASSERT(
                    typh_err == TYPH_SUCCESS,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Failed to delete old schedule");
        }

        typh_err = Build_Schedule(phase_id);
        TYPH_ASSERT(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Failed to build new schedule");

        phase->is_built = true;
    }

    // Commit needs doing once for pure Phases (on first pass)
    if (!phase->is_committed) {
        typh_err = Commit_Phase(phase);
        TYPH_ASSERT(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Failed to commit phase");

        phase->is_committed = true;
    }

    Schedule *schedule = phase->schedule;
    TYPH_ASSERT(
            schedule != nullptr,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Schedule missing after Build_Schedule");

    int const num_send = schedule->num_send;
    int const num_recv = schedule->num_recv;
    MPI_Request *recv_requests = schedule->recv_requests;
    int const mpi_tag = phase_id;

    // Post Receives
    if (num_recv > 0) {
        std::fill(recv_requests, recv_requests + num_recv, MPI_REQUEST_NULL);

        for (int i = 0; i < num_recv; i++) {
            if (schedule->mpi_recv_types[i] != MPI_DATATYPE_NULL) {
                int mpi_err =
                    MPI_Irecv(
                        MPI_BOTTOM,
                        1,
                        schedule->mpi_recv_types[i],
                        schedule->recv_procs[i],
                        mpi_tag,
                        mpi_runtime->comm,
                        &recv_requests[i]);
                TYPH_ASSERT(
                        mpi_err == MPI_SUCCESS,
                        ERR_MPI,
                        TYPH_ERR_MPI,
                        "MPI_Irecv failed");
            }
        }
    }

    // Perform Sends
    if (num_send > 0) {
        for (int i = 0; i < num_send; i++) {
            if (schedule->mpi_send_types[i] != MPI_DATATYPE_NULL) {
                int mpi_err =
                    MPI_Send(
                        MPI_BOTTOM,
                        1,
                        schedule->mpi_send_types[i],
                        schedule->send_procs[i],
                        mpi_tag,
                        mpi_runtime->comm);
                TYPH_ASSERT(
                        mpi_err == MPI_SUCCESS,
                        ERR_MPI,
                        TYPH_ERR_MPI,
                        "MPI_Send failed");
            }
        }
    }

    return TYPH_SUCCESS;
}



/**
 * \param [in] phase_id     the phase to finalise the exchange for
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
Finish_Exchange(int phase_id)
{
    Phase *phase;
    int irc = TYPH_REGISTRY->Get_Phase(phase_id, &phase);
    TYPH_ASSERT(
            irc == TYPH_SUCCESS,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Phase " + std::to_string(phase_id) + " not found");

    Schedule *schedule = phase->schedule;
    TYPH_ASSERT(
            schedule != nullptr,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Schedule absent for phase");

    int const num_recv = schedule->num_recv;
    if (num_recv > 0) {
        MPI_Request *recv_requests = schedule->recv_requests;
        TYPH_ASSERT(
                recv_requests != nullptr,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Receive requests absent for schedule");

        int mpi_err = MPI_Waitall(num_recv, recv_requests, MPI_STATUS_IGNORE);
        TYPH_ASSERT(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Waitall failed");
    }

    return irc;
}

} // namespace
} // namespace _TYPH_Internal

/**
 * \param [in] phase_id     the phase to exchange
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Exchange(int phase_id)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    int typh_err = Start_Exchange(phase_id);
    TYPH_ASSERT_RET(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Start_Exchange failed");

    return Finish_Exchange(phase_id);
}



/**
 * \param [in] phase_id     the phase to start the exchange for
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Start_Exchange(int phase_id)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    return Start_Exchange(phase_id);
}



/**
 * \param [in] phase_id     the phase to finalise the exchange for
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Finish_Exchange(int phase_id)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_CORE != nullptr && TYPH_CORE->Is_Initialised()),
            ERR_USER,
            TYPH_ERR_UNINITIALISED,
            "Not initialised");

    return Finish_Exchange(phase_id);
}
