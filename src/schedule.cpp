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
#include <iostream>
#include <numeric>
#include <cassert>
#include <memory>

#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "core.h"
#include "register.h"
#include "decomposition.h"
#include "keys.h"
#include "schedule.h"



namespace _TYPH_Internal {

int Schedule::max_parts = 0;
int Schedule::max_procs = 0;

namespace {

/**
 * \addtogroup internal
 *
 * @{
 */

/**
 * \brief   Check if there is any data to send/recv within the given key set and
 *          layers.
 *
 * \param [in]  key_set
 * \param [in]  sendrecv    send or recv keys?
 * \param [in]  min_layer   minimum ghost layer
 * \param [in]  max_layer   maximum ghost layer
 * \param [out] any         is there any data?
 *
 * \returns success or failure
 */
int
Any_Sendrecv(Key_Set const *key_set, TYPH_Sendrecv sendrecv, int min_layer,
        int max_layer, bool *any)
{
    *any = false;

    if (key_set == nullptr) {
        return TYPH_FAIL;
    }

    Key const *head = nullptr;
    switch (sendrecv) {
    case TYPH_SENDRECV_SEND: head = key_set->send_keys; break;
    case TYPH_SENDRECV_RECV: head = key_set->recv_keys; break;
    default:
        TYPH_ASSERT_RET(
                false,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Invalid sendrecv");
    }

    if (head != nullptr) {
        Key const *key = head;
        while (key != nullptr) {
            int const layer = key->layer;
            if ((layer >= min_layer) && (layer <= max_layer) &&
                    key->list_len > 0) {
                *any = true;
                break;
            }

            key = key->next;
        }
    }

    return TYPH_SUCCESS;
}



/**
 * \brief   Accumulate send/recv counts per processor from key set, between
 *          layers.
 *
 * \param [in]  key_set
 * \param [in]  sendrecv    send or recv keys?
 * \param [in]  min_layer   minimum ghost layer
 * \param [in]  max_layer   maximum ghost layer
 * \param [out] counts      per proc send/recv counts
 *
 * \return  success or failure
 */
int
Count_Sendrecv(Key_Set const *key_set, TYPH_Sendrecv sendrecv, int min_layer,
        int max_layer, int *counts)
{
    if (key_set == nullptr) {
        return TYPH_FAIL;
    }

    Key const *head = nullptr;
    switch (sendrecv) {
    case TYPH_SENDRECV_SEND: head = key_set->send_keys; break;
    case TYPH_SENDRECV_RECV: head = key_set->recv_keys; break;
    default:
        TYPH_ASSERT_RET(
                false,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Invalid sendrecv");
    }

    if (head != nullptr) {
        Key const *key = head;
        while (key != nullptr) {
            int const layer = key->layer;
            if ((layer >= min_layer) && (layer <= max_layer)) {
                counts[key->proc]++;
            }

            key = key->next;
        }
    }

    return TYPH_SUCCESS;
}



int
Is_Quant_Address_Set(int quant_id, bool *set)
{
    *set = false;

    Quant *quant = nullptr;
    int typh_err = TYPH_REGISTRY->Get_Quant(quant_id, &quant);
    if (typh_err != TYPH_SUCCESS) return TYPH_FAIL;

    TYPH_ASSERT_RET(
            quant != nullptr,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Quant null");

    *set = quant->quant_address != TYPH_NULL_ADDR;
    return TYPH_SUCCESS;
}



int
Populate_Quant_Info(Quant_Info *quant_info)
{
    Quant *quant = nullptr;
    int typh_err = TYPH_REGISTRY->Get_Quant(quant_info->quant_id, &quant);
    TYPH_ASSERT_RET(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Quant " + std::to_string(quant_info->quant_id) + " not found");

    // Check Ghost range is valid for sending Quant
    Key_Set const *key_set = quant_info->key_set;
    int const min_layer = quant_info->layer_min;
    int const max_layer = quant_info->layer_max;

    TYPH_ASSERT_RET(
            (key_set->layer_min <= min_layer) ||
                (key_set->layer_max >= max_layer),
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Key set has insufficient ghost layers for quant");

    // Set dimension for this quant
    int quant_size = 1;
    if (quant->rank > 1) {
        int *dims = new int[quant->rank - 1];
        for (int i = 0, j = 0; i < quant->rank; i++) {
            if (quant->dims[i] != TYPH_MESH_DIM) {
                dims[j++] = quant->dims[i];
            }
        }

        int const dim_prod = std::accumulate(dims, dims + quant->rank - 1,
                1, std::multiplies<int>());
        delete[] dims;
        quant_size = dim_prod / key_set->stride;
    }

    //
    // For Quants with >1D where index is not ordered (for example) (nel,3)
    // we set the quantsize to 1 and define the quant's stride as the mesh
    // dimension value.
    //
    // e.g.
    // arr(nel,3) - "default" arrangement, stride=1,  quantsize=3,nrepeat=1
    // arr(3,nel) - "reversed"             stride=nel,quantsize=1,nrepeat=1
    //
    // As we use the mesh based array bounds the stride will also be worked
    // out correctly where the mesh dimension is (for example)
    // arr(3,0:nel+2)
    //
    // Will also work for 3D arrays with "wrong" order, e.g. arr(3,nel,2)
    // quant->mesh_dim gives the dimension that is mesh based (1 in this
    // case)
    //
    if (quant->mesh_dim != 0) {
        quant_info->nrepeat = quant_size;
        quant_info->quant_size = 1;

        // Use upper/lower bound of array to work out stride
        quant_info->stride =
            quant->upper_bound[quant->mesh_dim] -
            quant->lower_bound[quant->mesh_dim];
        // e.g. arr(1:5)=>stride=5
        //      arr(0:5)=>stride=6

    } else {
        // nothing special needs doing for the "default" ordering
        quant_info->nrepeat = 1;
        quant_info->quant_size = quant_size;
        quant_info->stride = 1;
    }

    // this may be NULL type. Will work out later
    quant_info->old_mpi_type = &quant->mpi_datatype;
    return TYPH_SUCCESS;
}



/**
 * \brief   TODO
 *
 * \param [in]  quant_info  phase quant info
 * \param [in]  proc        target proc
 * \param [in]  sendrecv    send or recv keys?
 * \param [out] parts       array segment (preallocated) to store parts
 */
int
Build_Schedule_Parts(Quant_Info const *quant_info, int proc,
        TYPH_Sendrecv sendrecv, Schedule_Part *parts)
{
    int idx = 0;
    while (true) {
        if (quant_info == nullptr) {
            break;
        }

        // Get quant associated with quant info
        Quant *quant;
        int typh_err = TYPH_REGISTRY->Get_Quant(
                quant_info->quant_id,
                &quant);
        TYPH_ASSERT_RET(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Get_Quant failed");

        // Don't set schedule for this quant if it has a NULL address
        if (quant->quant_address == TYPH_NULL_ADDR) {
            quant_info = quant_info->next;
            continue;
        }

        int const min_layer = quant_info->layer_min;
        int const max_layer = quant_info->layer_max;
        for (int ilayer = min_layer; ilayer <= max_layer; ilayer++) {

            // Get key for this proc and layer
            int key_id;
            Key *key;
            typh_err = Get_Key(
                    quant_info->key_set,
                    sendrecv,
                    proc,
                    ilayer,
                    key_id,
                    &key);

            if (typh_err == TYPH_SUCCESS) {
                TYPH_ASSERT_RET(
                        key != nullptr,
                        ERR_INT,
                        TYPH_ERR_INTERNAL,
                        "Key null");

                parts[idx].key          = key;
                parts[idx].quant_size   = quant_info->quant_size;
                parts[idx].old_mpi_type = quant_info->old_mpi_type;
                parts[idx].address      = &quant->quant_address;
                parts[idx].nrepeat      = quant_info->nrepeat;
                parts[idx].stride       = quant_info->stride;
                parts[idx].new_mpi_type = MPI_DATATYPE_NULL;

                idx++;
            }
        }

        quant_info = quant_info->next;

    } // while true

    return TYPH_SUCCESS;
}



/**
 * \brief   TODO
 */
int
Commit_Schedule_Part(Schedule_Part *schedule_part)
{
    if (schedule_part->new_mpi_type != MPI_DATATYPE_NULL) {
        int mpi_err = MPI_Type_free(&schedule_part->new_mpi_type);
        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Type_free failed");
    }

    // If the key contains data...
    if (schedule_part->key->list_len > 0) {
        int const nlist =
            schedule_part->key->list_len * schedule_part->nrepeat;

        // For 2D/3D arrays with the element/node index at the start, we
        // need to create multiples of the displacement index. e.g.
        // (nel,3) arrays will effectively be three copies of a nel
        // arrays displacements whereas a (3,nel) array will have nel
        // sets of blocks of 3

        int *plist;
        if (schedule_part->nrepeat > 1) {
            plist = new int[nlist];

            int list = 0;
            for (int i = 0; i < schedule_part->nrepeat; i++) {
                for (int j = 0; j < schedule_part->key->list_len; j++) {
                    plist[list] =
                        i * schedule_part->stride +
                        schedule_part->key->list[j];
                    list++;
                }
            }

            TYPH_ASSERT_RET(
                    list == nlist,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "plist count mismatch");

        } else {
            plist = new int[nlist];
            std::copy(
                    schedule_part->key->list,
                    schedule_part->key->list + nlist,
                    plist);
        }

        if (schedule_part->key->block_lens != nullptr) {
            int mpi_err = MPI_Type_indexed(
                    nlist,
                    schedule_part->key->block_lens,
                    plist,
                    *schedule_part->old_mpi_type,
                    &schedule_part->new_mpi_type);

            TYPH_ASSERT_RET(
                    mpi_err == MPI_SUCCESS,
                    ERR_MPI,
                    TYPH_ERR_MPI,
                    "MPI_Type_indexed failed");

        } else {
            std::transform(plist, plist + nlist, plist,
                    [schedule_part](int v) {
                        return v * schedule_part->quant_size;
                    });

            int mpi_err = MPI_Type_create_indexed_block(
                    nlist,
                    schedule_part->quant_size,
                    plist,
                    *schedule_part->old_mpi_type,
                    &schedule_part->new_mpi_type);

            TYPH_ASSERT_RET(
                    mpi_err == MPI_SUCCESS,
                    ERR_MPI,
                    TYPH_ERR_MPI,
                    "MPI_Type_create_indexed_block failed");
        }

        delete[] plist;
    } // if (schedule_part->key->list_len > 0)

    return TYPH_SUCCESS;
}

/** @} */

} // namespace

/**
 * \param [in] phase_id     Phase ID we want to build the schedule for
 * \param [in] num_layers   Number of ghost layers
 *
 * \returns success or failure
 */
int
Build_Schedule(int phase_id, int num_layers __attribute__((unused)))
{
    int typh_err;
    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    // Get phase, and the phase's keyset to get a pointer to the MPI info
    Phase *phase;
    typh_err = TYPH_REGISTRY->Get_Phase(phase_id, &phase);
    TYPH_ASSERT_RET(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Phase " + std::to_string(phase_id) + " not found");

    Key_Set *key_set;
    typh_err = Get_Key_Set(phase->key_set_id, &key_set);
    TYPH_ASSERT_RET(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Key set " + std::to_string(phase->key_set_id) + " not found");

    Partition_Info *partition_info;
    typh_err = Get_Partition_Info(key_set->partition_id, &partition_info);
    TYPH_ASSERT_RET(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Partition info " + std::to_string(key_set->partition_id) +
                " not found");

    // These will store the number of keys that will be sent or received to/from
    // each proc (keysets may have keys for layers we don't care about so we'll
    // ignore the keys outside our range)
    std::unique_ptr<int[]> send_counts(new int[mp->size]);
    TYPH_ALLOC_CHECK(send_counts.get(), "send_counts");
    std::fill(send_counts.get(), send_counts.get() + mp->size, 0);

    std::unique_ptr<int[]> recv_counts(new int[mp->size]);
    TYPH_ALLOC_CHECK(recv_counts.get(), "recv_counts");
    std::fill(recv_counts.get(), recv_counts.get() + mp->size, 0);

    // Loop over all quants in this phase and get quant info for each
    Quant_Info *quant_info = phase->quant_info;
    for (int k = 0; k < phase->num_quants; k++) {
        TYPH_ASSERT_RET(
                quant_info != nullptr,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Fewer quants than expected");

        // Get Key_Set for this Quant in this Phase. Each Quant could have
        // unique Key_Set.
        typh_err = Get_Key_Set(quant_info->key_set_id, &key_set);
        TYPH_ASSERT_RET(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Key set " + std::to_string(quant_info->key_set_id) +
                    " not found");

        quant_info->key_set = key_set;

        // Ghost layers are set to the value stored in the quant (which default
        // to the phase values)
        int const min_layer = quant_info->layer_min;
        int const max_layer = quant_info->layer_max;

        // Send and receive Quant can be different. One or other could be
        // unallocated. We need to work out if a) anything is being sent or
        // received and b) is the array allocated

        // Part a) Anything being sent/received?
        // We don't need the exact number of sends/receives, just if it is > 0
        bool any_send, any_recv;
        typh_err  = Any_Sendrecv(key_set, TYPH_SENDRECV_SEND, min_layer,
                max_layer, &any_send);
        typh_err |= Any_Sendrecv(key_set, TYPH_SENDRECV_RECV, min_layer,
                max_layer, &any_recv);
        TYPH_ASSERT_RET(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Any_Sendrecv failed");

        // Now part b) Is the array associated?
        // Note: the schedule builder will ignore unassociated Quants that
        // should be sending/receiving data. Assumption is that this Quant is
        // not involved.  This can occur for Phases where the Quant can change
        // each time This can lead to errors though where the user forgot to set
        // the Quant address for a valid Quant Unassociated Quants that are not
        // sending/receiving will be ignored later on
        bool send_quant_allocd, recv_quant_allocd;
        typh_err  = Is_Quant_Address_Set(quant_info->quant_id,
                &send_quant_allocd);
        typh_err |= Is_Quant_Address_Set(quant_info->recv_quant_id,
                &recv_quant_allocd);
        TYPH_ASSERT_RET(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Is_Quant_Address_Set failed");

        if ((!send_quant_allocd && any_send) ||
                (!recv_quant_allocd && any_recv)) {
            quant_info = quant_info->next;
            continue;
        }

        // Determine quant info size, stride, nrepeat etc.
        typh_err = Populate_Quant_Info(quant_info);
        TYPH_ASSERT_RET(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Populate_Quant_Info failed");

        // Work out exact number of send and receive keys to/from each
        // neighbouring processor
        if (send_quant_allocd) {
            typh_err = Count_Sendrecv(key_set, TYPH_SENDRECV_SEND,
                    min_layer, max_layer, send_counts.get());
            TYPH_ASSERT_RET(
                    typh_err == TYPH_SUCCESS,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Count_Sendrecv failed");
        }

        if (recv_quant_allocd) {
            typh_err = Count_Sendrecv(key_set, TYPH_SENDRECV_RECV, min_layer,
                    max_layer, recv_counts.get());
            TYPH_ASSERT_RET(
                    typh_err == TYPH_SUCCESS,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Count_Sendrecv failed");
        }

        // Next Quant in this phase
        quant_info = quant_info->next;

    } // for 0 <= k < phase->num_quants

    if (quant_info != nullptr) {
        TYPH_WARNING("More quants than expected");
    }

    // Allocate new schedule
    Schedule *schedule = phase->schedule;
    TYPH_ASSERT_RET(
            schedule == nullptr,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Schedule already exists for phase");

    schedule = new Schedule();
    TYPH_ALLOC_CHECK(schedule, "schedule");
    phase->schedule = schedule;

    // Store the number of procs we are sending/receiving to/from
    schedule->num_send = 0;
    schedule->num_recv = 0;
    for (int proc = 0; proc < mp->size; proc++) {
        if (send_counts[proc] > 0) schedule->num_send++;
        if (recv_counts[proc] > 0) schedule->num_recv++;
    }

    // Allocate schedule parts
    int const total_count =
        std::accumulate(send_counts.get(), send_counts.get() + mp->size, 0) +
        std::accumulate(recv_counts.get(), recv_counts.get() + mp->size, 0);
    if (total_count > 0) {
        schedule->parts = new Schedule_Part[total_count];
        TYPH_ALLOC_CHECK(schedule->parts, "schedule->parts");
    }

    // Set up send data start positions and displacements. Each proc's messages
    // will be broken into nparts where nparts = no. send keys to that proc.
    int pos = 0;
    int num_send = schedule->num_send;
    if (num_send > 0) {
        schedule->send_procs = new int[num_send];
        schedule->mpi_send_types = new MPI_Datatype[num_send];
        schedule->send_requests = new MPI_Request[num_send];
        schedule->send_nparts = new int[num_send];
        schedule->send_start = new int[num_send];

        TYPH_ALLOC_CHECK(schedule->send_procs, "schedule->send_procs");
        TYPH_ALLOC_CHECK(schedule->mpi_send_types, "schedule->mpi_send_types");
        TYPH_ALLOC_CHECK(schedule->send_requests, "schedule->send_requests");
        TYPH_ALLOC_CHECK(schedule->send_nparts, "schedule->send_nparts");
        TYPH_ALLOC_CHECK(schedule->send_start, "schedule->send_start");

        int idx = 0;
        for (int iproc = 0; iproc < mp->size; iproc++) {
            if (send_counts[iproc] > 0) {
                schedule->send_procs[idx]     = iproc;
                schedule->mpi_send_types[idx] = MPI_DATATYPE_NULL;
                schedule->send_nparts[idx]    = send_counts[iproc];
                schedule->send_start[idx]     = pos;

                idx++;
                pos += send_counts[iproc];
            }
        }
        TYPH_ASSERT_RET(
                idx == num_send,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "num_send value mismatch");

        // Loop over quants and then ghost layers within each quant. Add each
        // valid send key as we find it.
        for (int isend = 0; isend < schedule->num_send; isend++) {
            int const proc  = schedule->send_procs[isend];
            int       start = schedule->send_start[isend];

            Build_Schedule_Parts(phase->quant_info, proc,
                    TYPH_SENDRECV_SEND, &schedule->parts[start]);
        }
    } // num_send > 0

    // Set up receive data start positions and displacements
    int num_recv = schedule->num_recv;
    if (num_recv > 0) {
        schedule->recv_procs = new int[num_recv];
        schedule->mpi_recv_types = new MPI_Datatype[num_recv];
        schedule->recv_requests = new MPI_Request[num_recv];
        schedule->recv_nparts = new int[num_recv];
        schedule->recv_start = new int[num_recv];

        TYPH_ALLOC_CHECK(schedule->recv_procs, "schedule->recv_procs");
        TYPH_ALLOC_CHECK(schedule->mpi_recv_types, "schedule->mpi_recv_types");
        TYPH_ALLOC_CHECK(schedule->recv_requests, "schedule->recv_requests");
        TYPH_ALLOC_CHECK(schedule->recv_nparts, "schedule->recv_nparts");
        TYPH_ALLOC_CHECK(schedule->recv_start, "schedule->recv_start");

        int idx = 0;
        for (int iproc = 0; iproc < mp->size; iproc++) {
            if (recv_counts[iproc] > 0) {
                schedule->recv_procs[idx]     = iproc;
                schedule->mpi_recv_types[idx] = MPI_DATATYPE_NULL;
                schedule->recv_nparts[idx]    = recv_counts[iproc];
                schedule->recv_start[idx]     = pos;

                idx++;
                pos += recv_counts[iproc];
            }
        }
        TYPH_ASSERT_RET(
                idx == num_recv,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "num_recv value mismatch");

        // Loop over Quants and then ghost layers within each Quant. Add each
        // valid recv key as we find it.
        for (int irecv = 0; irecv < schedule->num_recv; irecv++) {
            int const proc  = schedule->recv_procs[irecv];
            int       start = schedule->recv_start[irecv];

            Build_Schedule_Parts(phase->quant_info, proc,
                    TYPH_SENDRECV_RECV, &schedule->parts[start]);
        }
    } // num_recv > 0

    Schedule::max_procs = std::max(Schedule::max_procs,
            std::max(schedule->num_send, schedule->num_recv));

    int const max_send_counts =
        *std::max_element(send_counts.get(), send_counts.get() + mp->size);
    int const max_recv_counts =
        *std::max_element(recv_counts.get(), recv_counts.get() + mp->size);
    Schedule::max_parts = std::max(Schedule::max_parts,
            std::max(max_send_counts, max_recv_counts));

    return TYPH_SUCCESS;
}



int
Delete_Schedule(int phase_id)
{
    int typh_err = TYPH_SUCCESS;
    int  mpi_err = MPI_SUCCESS;

    Phase *phase;
    typh_err = TYPH_REGISTRY->Get_Phase(phase_id, &phase);
    TYPH_ASSERT_RET(
            typh_err == TYPH_SUCCESS,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Failed to find phase from phase_id");

    Schedule *schedule = phase->schedule;
    if (schedule == nullptr) {
        return TYPH_SUCCESS;
    }

    // Free MPI types
    if (schedule->mpi_send_types != nullptr) {
        for (int i = 0; i < schedule->num_send; i++) {
            if (schedule->mpi_send_types[i] != MPI_DATATYPE_NULL) {
                mpi_err |= MPI_Type_free(&schedule->mpi_send_types[i]);
            }
        }
    }

    if (schedule->mpi_recv_types != nullptr) {
        for (int i = 0; i < schedule->num_recv; i++) {
            if (schedule->mpi_recv_types[i] != MPI_DATATYPE_NULL) {
                mpi_err |= MPI_Type_free(&schedule->mpi_recv_types[i]);
            }
        }
    }

    if (schedule->parts != nullptr) {
        for (int i = 0; i < schedule->num_send + schedule->num_recv; i++) {
            Schedule_Part *schedule_part = &schedule->parts[i];
            if (schedule_part != nullptr) {
                if (schedule_part->new_mpi_type != MPI_DATATYPE_NULL) {
                    mpi_err |= MPI_Type_free(&schedule_part->new_mpi_type);
                }
            }

            schedule_part->key          = nullptr;
            schedule_part->address      = nullptr;
            schedule_part->old_mpi_type = nullptr;
        }
    }

    TYPH_ASSERT_RET(
            mpi_err == MPI_SUCCESS,
            ERR_MPI,
            TYPH_ERR_MPI,
            "MPI_Type_free failed");

    // Deallocate everything inside the schedule
    #define SAFE_DELETE(ptr) { \
        if ((ptr) != nullptr) { \
            delete[] (ptr); \
            (ptr) = nullptr; \
        }}

    SAFE_DELETE(schedule->send_procs);
    SAFE_DELETE(schedule->mpi_send_types);
    SAFE_DELETE(schedule->send_requests);
    SAFE_DELETE(schedule->send_nparts);
    SAFE_DELETE(schedule->send_start);

    SAFE_DELETE(schedule->recv_procs);
    SAFE_DELETE(schedule->mpi_recv_types);
    SAFE_DELETE(schedule->recv_requests);
    SAFE_DELETE(schedule->recv_nparts);
    SAFE_DELETE(schedule->recv_start);

    SAFE_DELETE(schedule->parts);

    #undef SAFE_DELETE

    // Then deallocate the schedule itself
    delete schedule;
    phase->schedule = nullptr;

    return TYPH_SUCCESS;
}



/**
 * \param [inout]   phase   phase to commit
 *
 * \returns success or failure
 */
int
Commit_Phase(Phase *phase)
{
    std::unique_ptr<MPI_Datatype[]> mpi_types(
            new MPI_Datatype[Schedule::max_parts]);
    TYPH_ALLOC_CHECK(mpi_types.get(), "mpi_types");
    std::fill(mpi_types.get(), mpi_types.get() + Schedule::max_parts,
            MPI_DATATYPE_NULL);

    std::unique_ptr<int[]> mpi_block_lens(new int[Schedule::max_parts]);
    TYPH_ALLOC_CHECK(mpi_block_lens.get(), "mpi_block_lens");
    std::fill(mpi_block_lens.get(), mpi_block_lens.get() + Schedule::max_parts,
            0);

    std::unique_ptr<MPI_Aint[]> mpi_addresses(
            new MPI_Aint[Schedule::max_parts]);
    TYPH_ALLOC_CHECK(mpi_addresses.get(), "mpi_addresses");
    std::fill(mpi_addresses.get(), mpi_addresses.get() + Schedule::max_parts,
            TYPH_NULL_ADDR);

    Schedule *schedule = phase->schedule;
    TYPH_ASSERT_RET(
            schedule != nullptr,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Schedule not built");

    // *** Sends ***
    for (int isend = 0; isend < schedule->num_send; isend++) {

        // Ensure the type has been freed up
        if (schedule->mpi_send_types[isend] != MPI_DATATYPE_NULL) {
            int mpi_err = MPI_Type_free(&schedule->mpi_send_types[isend]);
            TYPH_ASSERT_RET(
                    mpi_err == MPI_SUCCESS,
                    ERR_MPI,
                    TYPH_ERR_MPI,
                    "MPI_Type_free failed");
        }

        // Build loop over destination procs
        int count = 0;
        int const send_end =
            schedule->send_start[isend] + schedule->send_nparts[isend];
        for (int iproc = schedule->send_start[isend]; iproc < send_end;
                iproc++) {

            Schedule_Part *schedule_part = &schedule->parts[iproc];
            int typh_err = Commit_Schedule_Part(schedule_part);
            TYPH_ASSERT_RET(
                    typh_err == TYPH_SUCCESS,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Commit_Schedule_Part failed");

            // Store information needed to create the final derived type
            if (schedule_part->key->list_len > 0) {
                mpi_types[count]      = schedule_part->new_mpi_type;
                mpi_block_lens[count] = 1;
                mpi_addresses[count]  = *schedule_part->address;
                count++;
            }
        }

        // Create derived type and commit
        if (count > 0) {
            int mpi_err = MPI_Type_create_struct(
                    count,
                    mpi_block_lens.get(),
                    mpi_addresses.get(),
                    mpi_types.get(),
                    &schedule->mpi_send_types[isend]);

            mpi_err |= MPI_Type_commit(&schedule->mpi_send_types[isend]);
            TYPH_ASSERT_RET(
                    mpi_err == MPI_SUCCESS,
                    ERR_MPI,
                    TYPH_ERR_MPI,
                    "Failed to create phase derived type");
        }
    } // for isend

    // *** Receives ***
    for (int irecv = 0; irecv < schedule->num_recv; irecv++) {

        // Ensure the type has been freed up
        if (schedule->mpi_recv_types[irecv] != MPI_DATATYPE_NULL) {
            int mpi_err = MPI_Type_free(&schedule->mpi_recv_types[irecv]);
            TYPH_ASSERT_RET(
                    mpi_err == MPI_SUCCESS,
                    ERR_MPI,
                    TYPH_ERR_MPI,
                    "MPI_Type_free failed");
        }

        // Build loop over destination procs
        int count = 0;
        int const recv_end =
            schedule->recv_start[irecv] + schedule->recv_nparts[irecv];
        for (int iproc = schedule->recv_start[irecv]; iproc < recv_end; iproc++) {

            Schedule_Part *schedule_part = &schedule->parts[iproc];
            int typh_err = Commit_Schedule_Part(schedule_part);
            TYPH_ASSERT_RET(
                    typh_err == TYPH_SUCCESS,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Commit_Schedule_Part failed");

            // Store information needed to create the final derived type
            if (schedule_part->key->list_len > 0) {
                mpi_types[count]      = schedule_part->new_mpi_type;
                mpi_block_lens[count] = 1;
                mpi_addresses[count]  = *schedule_part->address;
                count++;
            }
        }

        // Create derived type and commit
        if (count > 0) {
            int mpi_err = MPI_Type_create_struct(
                    count,
                    mpi_block_lens.get(),
                    mpi_addresses.get(),
                    mpi_types.get(),
                    &schedule->mpi_recv_types[irecv]);

            mpi_err |= MPI_Type_commit(&schedule->mpi_recv_types[irecv]);
            TYPH_ASSERT_RET(
                    mpi_err == MPI_SUCCESS,
                    ERR_MPI,
                    TYPH_ERR_MPI,
                    "Failed to create phase derived type");
        }
    } // for irecv

    return TYPH_SUCCESS;
} // Commit_Phase



int
Validate_Committed_Phase(Phase *phase)
{
    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    Schedule *schedule = phase->schedule;
    if (schedule == nullptr) {
        return TYPH_FAIL;
    }

    // Check that each send's size is consistent with the receiver's size
    std::unique_ptr<int[]> send_sizes(new int[schedule->num_send]);
    std::unique_ptr<int[]> recv_sizes(new int[schedule->num_recv]);

    std::fill(send_sizes.get(), send_sizes.get() + schedule->num_send, 0);
    std::fill(recv_sizes.get(), recv_sizes.get() + schedule->num_recv, 0);

    for (int isend = 0; isend < schedule->num_send; isend++) {
        int mpi_err = MPI_Type_size(
                schedule->mpi_send_types[isend], &send_sizes[isend]);
        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "Failed to measure MPI type size");
    }

    for (int irecv = 0; irecv < schedule->num_recv; irecv++) {
        int mpi_err = MPI_Type_size(
                schedule->mpi_recv_types[irecv], &recv_sizes[irecv]);
        TYPH_ASSERT_RET(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "Failed to measure MPI type size");
    }

    std::unique_ptr<int[]> sizes_from_sendr(new int[schedule->num_recv]);
    std::unique_ptr<int[]> sizes_from_recvr(new int[schedule->num_send]);

    std::unique_ptr<MPI_Request[]> recv_reqs(new MPI_Request[schedule->num_recv]);
    std::unique_ptr<MPI_Request[]> send_reqs(new MPI_Request[schedule->num_send]);

    std::fill(recv_reqs.get(), recv_reqs.get() + schedule->num_recv,
            MPI_REQUEST_NULL);
    std::fill(send_reqs.get(), send_reqs.get() + schedule->num_send,
            MPI_REQUEST_NULL);

    // ... Receive sent message sizes from sender
    for (int irecv = 0; irecv < schedule->num_recv; irecv++) {
        int mpi_err = MPI_Irecv(&sizes_from_sendr[irecv], 1, MPI_INTEGER,
                schedule->recv_procs[irecv], 0, mp->comm, &recv_reqs[irecv]);
        TYPH_ASSERT(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Irecv failed");
    }

    // ... Receive received message sizes from receiver
    for (int isend = 0; isend < schedule->num_send; isend++) {
        int mpi_err = MPI_Irecv(&sizes_from_recvr[isend], 1, MPI_INTEGER,
                schedule->send_procs[isend], 0, mp->comm, &send_reqs[isend]);
        TYPH_ASSERT(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Irecv failed");
    }

    // ... Send sent message sizes to receiver
    for (int isend = 0; isend < schedule->num_send; isend++) {
        int mpi_err = MPI_Send(&send_sizes[isend], 1, MPI_INTEGER,
                schedule->send_procs[isend], 0, mp->comm);
        TYPH_ASSERT(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Send failed");
    }

    // ... Send received message sizes to sender
    for (int irecv = 0; irecv < schedule->num_recv; irecv++) {
        int mpi_err = MPI_Send(&recv_sizes[irecv], 1, MPI_INTEGER,
                schedule->recv_procs[irecv], 0, mp->comm);
        TYPH_ASSERT(
                mpi_err == MPI_SUCCESS,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Send failed");
    }

    int mpi_err = MPI_Waitall(schedule->num_send, send_reqs.get(),
            MPI_STATUS_IGNORE);
    TYPH_ASSERT(
            mpi_err == MPI_SUCCESS,
            ERR_MPI,
            TYPH_ERR_MPI,
            "MPI_Waitall failed");

    mpi_err = MPI_Waitall(schedule->num_recv, recv_reqs.get(),
            MPI_STATUS_IGNORE);
    TYPH_ASSERT(
            mpi_err == MPI_SUCCESS,
            ERR_MPI,
            TYPH_ERR_MPI,
            "MPI_Waitall failed");

    // The size of the messages we send should match the size the receiver is
    // expecting
    for (int isend = 0; isend < schedule->num_send; isend++) {
        bool const correct = (send_sizes[isend] == sizes_from_recvr[isend]);
        TYPH_ASSERT_RET(
                correct,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Message sizes don't match between sender " +
                    std::to_string(mp->rank) + " and receiver " +
                    std::to_string(schedule->send_procs[isend]));
    }

    // The size of the messages we expect to receive should match the size the
    // sender is going to send
    for (int irecv = 0; irecv < schedule->num_recv; irecv++) {
        bool const correct = (recv_sizes[irecv] == sizes_from_sendr[irecv]);
        TYPH_ASSERT_RET(
                correct,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Message sizes don't match between receiver " +
                    std::to_string(mp->rank) + " and sender " +
                    std::to_string(schedule->recv_procs[irecv]));
    }

    return TYPH_SUCCESS;
}

} // namespace _TYPH_Internal
