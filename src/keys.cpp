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
#include <vector>
#include <cassert>
#include <memory>
#include <mpi.h>

#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "core.h"
#include "decomposition.h"
#include "keys.h"




namespace _TYPH_Internal {
namespace {

/**
 * \addtogroup internal
 *
 * @{
 */

/** Store created key sets. */
std::vector<Key_Set *> key_sets;

/**
 * \brief   Add a key set to the list.
 *
 * \param [out] key_set_id  the id of the added key set
 * \param [in]  key_set     the key set to add
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
Add_Key_Set(
        int &key_set_id,
        Key_Set *key_set)
{
    key_sets.push_back(key_set);
    key_set_id = key_sets.size() - 1;
    return TYPH_SUCCESS;
}



/**
 * \brief TODO
 */
int
New_Key_Set(int &key_set_id, TYPH_Centring centring, TYPH_Auxiliary aux,
        int stride, int layer_min, int layer_max, int partition_id,
        Key_Set *&key_set)
{
    using namespace _TYPH_Internal;

    key_set = new Key_Set();
    TYPH_ALLOC_CHECK(key_set, "key_set");

    // Add new keyset to keyset array
    int typh_err = Add_Key_Set(key_set_id, key_set);
    TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Add_Key_Set failed");

    key_set->centring     = centring;
    key_set->aux          = aux;
    key_set->stride       = stride;
    key_set->layer_min    = layer_min;
    key_set->layer_max    = layer_max;
    key_set->partition_id = partition_id;
    key_set->num_send     = 0;
    key_set->num_recv     = 0;
    key_set->send_keys    = nullptr;
    key_set->recv_keys    = nullptr;

    return TYPH_SUCCESS;
}



/**
 * \brief TODO
 */
int
Add_Key(Key *&head, int &num_keys, int proc, int layer, Key *&node)
{
    using namespace _TYPH_Internal;

    node = nullptr;

    // Keys are sorted by proc and then layer, so find the appropriate insertion
    // point
    Key *prev_key = nullptr;
    Key *cur_key = head;
    for (int i = 0; i < num_keys; i++) {
        TYPH_ASSERT(
                cur_key != nullptr,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Fewer keys than expected");

        if ((cur_key->proc > proc) ||
            (cur_key->proc == proc && cur_key->layer > layer)) {
            break;
        }

        prev_key = cur_key;
        cur_key = cur_key->next;
    }

    Key *new_key = new Key();
    TYPH_ALLOC_CHECK(new_key, "new_key");

    // Set next (nullptr if last item in the list)
    new_key->next = cur_key;

    // Set prev
    if (prev_key != nullptr) {
        prev_key->next = new_key;
        new_key->prev = prev_key;

    } else {
        head = new_key;
        new_key->prev = nullptr;
    }

    node = new_key;
    node->proc     = proc;
    node->layer    = layer;
    node->list_len = 0;
    node->parent   = nullptr;
    node->list     = nullptr;

    num_keys++;
    return TYPH_SUCCESS;
}



/**
 * \brief TODO
 */
int
Define_Key_Set(
        Key_Set *key_set,
        TYPH_Keytype key_type,
        int layer_min,
        int layer_max,
        TYPH_Centring centring __attribute__((unused)),
        int stride,
        int const *to_proc,
        int const *to_glob,
        int const *total_items,
        int const *connectivity)
{
    //centring unused??

    using namespace _TYPH_Internal;

    auto constexpr IXp = &Index_2D<2>;

    MPI_Status mpi_status;
    int mpi_err, typh_err;

    Key *key = nullptr;

    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    int proc = mp->rank;

    //integer(kind=INK) :: iStartindex  ! start index
    //integer(kind=INK) :: iLayer       ! current ghost layer
    //integer(kind=INK) :: iNsize       ! No. properties per communicated item (nd)
    //integer(kind=INK) :: iNitems      ! No. send items
    //integer(kind=INK) :: iNrecv       ! No. recv items
    //integer(kind=INK) :: iNcount      ! message count
    //integer(kind=INK) :: iNremain     ! remaining message count
    //integer(kind=INK) :: iSreq        ! send request counter

    int item_size;      // # properties per communicated item
    int num_send;       // # send items
    int num_recv;       // # recv items
    int num_msg;        // # messages
    int msg_remain;     // # messages remaining
    int send_req_count = 0; // Send request counter

    MPI_Request recv_req;
    int *send_item_list = nullptr;
    int *recv_item_list = nullptr;

    // If global map array is present ensure that this will be included in
    // send_item_list
    bool glob_present;
    if (to_glob != nullptr) {
        item_size = 2;
        glob_present = true;
    } else {
        item_size = 1;
        glob_present = false;
    }

    // # ghost items from each neighbour
    // 0..size-1
    std::unique_ptr<int[]> counts(new int[mp->size]);
    TYPH_ALLOC_CHECK(counts.get(), "counts");
    std::fill(counts.get(), counts.get() + mp->size, 0);

    // 0..size-1
    std::unique_ptr<int[]> starts(new int[mp->size]);
    TYPH_ALLOC_CHECK(starts.get(), "starts");
    std::fill(starts.get(), starts.get() + mp->size, 0);

    // 0..size-1
    std::unique_ptr<int[]> t_counts(new int[mp->size]);
    TYPH_ALLOC_CHECK(t_counts.get(), "t_counts");
    std::fill(t_counts.get(), t_counts.get() + mp->size, 0);

    // 1..size
    std::unique_ptr<MPI_Request[]> send_req(new MPI_Request[mp->size]);
    TYPH_ALLOC_CHECK(send_req.get(), "send_req");
    std::fill(send_req.get(), send_req.get() + mp->size, MPI_REQUEST_NULL);

    // Build key lists for each ghost cell in each ghost layer
    for (int layer = layer_min; layer <= layer_max; layer++) {

        #define IXsil(i, j) (Index_2D(i, j, item_size))

        std::fill(counts.get(), counts.get() + mp->size, 0);
        int const comm_tag = 1000 + layer;

        switch (key_type) {
        case TYPH_KEYTYPE_CELL:
        case TYPH_KEYTYPE_NODE:
            {
                num_send = 0;

                int start_index = 0;
                if (layer > 0) {
                    start_index = total_items[layer-1];
                }

                // First work out how many ghost cells/nodes belong to each
                // neighbour processor
                for (int i = start_index; i < total_items[layer]; i++) {
                    int const proc = to_proc[IXp(0, i)];
                    if (proc == mp->rank) {
                        continue;
                    }

                    num_send++;
                    counts[proc]++;
                }

                // Now build up a list of start positions for each procs
                // contributions in the send array that will be built shortly
                starts[0] = 0;
                for (int i = 1; i < mp->size; i++) {
                    starts[i] = starts[i-1] + counts[i-1];
                }

                // Build an array containing the displacements (localID-1) and
                // global index (if avail) of each ghost cell receiving data in
                // the exchange on this processor. iStarts tells us where each
                // neighbour proc's contribution begins. The relevent bits
                // passed to the neighbour procs later to let them work out
                // their sending cell displacements.

                // NOTE: this array can be used to work out a reverse send signal

                //allocate(iSendItemList(iNsize, iNitems), stat=iAllocStat)
                send_item_list = new int[item_size * num_send];
                TYPH_ALLOC_CHECK(send_item_list, "send_item_list");

                // Now loop through and fill the array
                // Could this be done without the need for the full second loop?
                for (int i = start_index; i < total_items[layer]; i++) {
                    int const proc = to_proc[IXp(0, i)];
                    if (proc == mp->rank) {
                        continue;
                    }

                    send_item_list[IXsil(0, starts[proc])] = i;
                    starts[proc]++; // Increases index in send_item_list
                }
            }
            break;

        case TYPH_KEYTYPE_CELL_CORNER:
            {
                TYPH_ASSERT(
                        connectivity != nullptr,
                        ERR_INT,
                        TYPH_ERR_INTERNAL,
                        "Connectivity required for cell corner key type");

                TYPH_ASSERT(
                        layer > 0,
                        ERR_INT,
                        TYPH_ERR_INTERNAL,
                        "Real layer invalid");

                Partition_Info *partition_info;
                typh_err = Get_Partition_Info(key_set->partition_id, &partition_info);
                TYPH_ASSERT(
                        typh_err == TYPH_SUCCESS,
                        ERR_INT,
                        TYPH_ERR_INTERNAL,
                        "Partition " + std::to_string(key_set->partition_id) +
                            " not found");

                // Cell corners are generated the same way with the added
                // complication that we need each corner's component

                num_send = 0;

                int start_index = total_items[layer-1];

                // First work out how many ghost cells/nodes belong to each
                // neighbour processor
                //int const j = partition_info->num_nd_total[layer-1];
                for (int i = start_index; i < total_items[layer]; i++) {
                    int const proc = to_proc[IXp(0, i)];
                    if (proc == mp->rank) {
                        continue;
                    }

                    for (int k = 0; k < partition_info->nodes_per_el; k++) {
                        //if (mConnectivity(kk, ii) <= jj) then
                        assert(false && "not yet implemented");
                        if (false) {
                            num_send++;
                            counts[proc]++;
                        }
                    }
                }

                // Now build up a list of start positions for each procs
                // contributions in the send array that will be built shortly
                starts[0] = 0;
                for (int i = 1; i < mp->size; i++) {
                    starts[i] = starts[i-1] + counts[i-1];
                }

                // Build an array containing the displacements (localID-1) and
                // global index (if avail) of each ghost cell receiving data in
                // the exchange on this processor. iStarts tells us where each
                // neighbour proc's contribution begins. The relevent bits
                // passed to the neighbour procs later to let them work out
                // their sending cell displacements.

                // NOTE: this array can be used to work out a reverse send signal

                send_item_list = new int[item_size * num_send];
                TYPH_ALLOC_CHECK(send_item_list, "send_item_list");

                // Now loop through and fill the array
                // Could this be done without the need for the full second loop?
                for (int i = start_index; i < total_items[layer]; i++) {
                    int const proc = to_proc[IXp(0, i)];
                    if (proc == mp->rank) {
                        continue;
                    }

                    for (int k = 0; k < partition_info->nodes_per_el; k++) {
                        //if (mConnectivity(kk, ii) <= jj) then
                        assert(false && "not yet implemented");
                        if (false) {
                            send_item_list[IXsil(0, starts[proc])] = i * partition_info->nodes_per_el + k;
                            starts[proc]++; // Increases index in send_item_list
                        }
                    }
                }
            }
            break;

        default:
            TYPH_ASSERT(false, ERR_INT, TYPH_ERR_INTERNAL, "Invalid key type");
        } // switch key_type

        // Reset starts array (was mangled in above loop constructing
        // send_item_list)
        starts[0] = 0;
        for (int i = 1; i < mp->size; i++) {
            starts[i] = starts[i-1] + counts[i-1];
        }

        // Ensure all processors have the total number of items to send

        // Must ensure this call is maintained if the collectives change
        typh_err = TYPH_Reduce_i(
                counts.get(),
                &mp->size,
                1,
                t_counts.get(),
                TYPH_OP_SUM);
        TYPH_ASSERT(typh_err == TYPH_SUCCESS, ERR_INT, TYPH_ERR_INTERNAL,
                "TYPH_Reduce_i failed");
        num_recv = t_counts[mp->rank];

        // Fill send_item_list with global IDs if available
        if (glob_present) {
            for (int i = 0; i < num_send; i++) {
                int const jj = send_item_list[IXsil(0, i)] / stride;
                send_item_list[IXsil(1, i)] = to_glob[jj];
            }
        }

        // Now have the information required to build the receive MPI types for
        // each processor. Time to build the receive keys and work out the send
        // information on the neighbour processors. Each processor may be
        // sending to more than one neighbour so need to build up list of all
        // communications
        send_req_count = 0;

        // Build receive key lists for the current processor from each
        // neighbour. Then convert send_item_list array to store local
        // displacement on neighbour proc. Then send out sections of the
        // send_item_list to the relevent neighbour proc.

        // First post a recv to prevent unexpected messages later
        recv_item_list = new int[item_size * num_recv];
        TYPH_ALLOC_CHECK(recv_item_list, "recv_item_list");
        #define IXril(i, j) (Index_2D(i, j, item_size))

        if (num_recv > 0) {
            mpi_err = MPI_Irecv(recv_item_list, item_size * num_recv,
                    MPI_INT, MPI_ANY_SOURCE, comm_tag, mp->comm, &recv_req);
            TYPH_ASSERT(mpi_err == MPI_SUCCESS, ERR_MPI, TYPH_ERR_MPI,
                    "MPI_Irecv failed");
        }

        for (int proc = 0; proc < mp->size; proc++) {
            if (counts[proc] > 0) {
                typh_err = Add_Key(
                        key_set->recv_keys,
                        key_set->num_recv,
                        proc,
                        layer,
                        key);
                TYPH_ASSERT(
                        typh_err == TYPH_SUCCESS,
                        ERR_INT,
                        TYPH_ERR_INTERNAL,
                        "Add_Key failed");

                key->list_len = counts[proc];

                key->list = new int[key->list_len];
                TYPH_ALLOC_CHECK(key->list, "key->list");

                int ii = 0;
                for (int jj = starts[proc]; jj < starts[proc] + counts[proc];
                        jj++) {

                    // Build receive key lists for the current processor from
                    // each neighbour.
                    key->list[ii] = send_item_list[IXsil(0, jj)];

                    int const mm = send_item_list[IXsil(0, jj)] / stride;
                    int const kk = send_item_list[IXsil(0, jj)] - (stride * mm);

                    // Convert send_item_list array to store local displacement
                    // on neighbour proc
                    send_item_list[IXsil(0, jj)] = stride * to_proc[IXp(1, mm)] + kk;
                    ii++;
                }

                // Send out sections of the send_item_list to the relevent
                // neighbour proc
                mpi_err = MPI_Isend(
                        &send_item_list[IXsil(0, starts[proc])],
                        item_size * counts[proc],
                        MPI_INT,
                        proc,
                        comm_tag,
                        mp->comm,
                        &send_req[send_req_count++]);
                TYPH_ASSERT(mpi_err == MPI_SUCCESS, ERR_MPI, TYPH_ERR_MPI,
                        "MPI_Isend failed");

            } // if (counts[proc] > 0)
        } // foreach proc

        // Receive sections of the send_item_list from neighbour procs
        msg_remain = num_recv;
        while (true) {
            if (msg_remain == 0) break;

            mpi_err = MPI_Wait(&recv_req, &mpi_status);
            TYPH_ASSERT(
                    mpi_err == MPI_SUCCESS,
                    ERR_MPI,
                    TYPH_ERR_MPI,
                    "MPI_Wait failed");

            mpi_err = MPI_Get_count(&mpi_status, MPI_INT, &num_msg);
            TYPH_ASSERT(
                    mpi_err == MPI_SUCCESS,
                    ERR_MPI,
                    TYPH_ERR_MPI,
                    "MPI_Get_count failed");

            num_msg /= item_size;
            msg_remain = msg_remain - num_msg;
            proc = mpi_status.MPI_SOURCE;

            typh_err = Add_Key(
                    key_set->send_keys,
                    key_set->num_send,
                    proc,
                    layer,
                    key);
            TYPH_ASSERT(
                    typh_err == TYPH_SUCCESS,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Add_Key failed");

            key->list_len = num_msg;

            key->list = new int[key->list_len];
            TYPH_ALLOC_CHECK(key->list, "key->list");

            // Build send key lists from the current processor to each neighbour
            for (int i = 0; i < num_msg; i++) {
                key->list[i] = recv_item_list[IXril(0, i)];
            }

            // Sanity check. Do the recv_item_list global IDs match the main
            // global list ones?
            if (glob_present) {
                for (int ii = 0; ii < num_msg; ii++) {
                    int const jj = recv_item_list[IXril(0, ii)] / stride;
                    if (to_glob[jj] != recv_item_list[IXril(1, ii)]) {
                        TYPH_ASSERT(false, ERR_INT, TYPH_ERR_INTERNAL,
                                "Global IDs do not match between processors");
                    }
                }
            }

            if (msg_remain > 0) {
                mpi_err = MPI_Irecv(
                        recv_item_list,
                        item_size * num_recv,
                        MPI_INT,
                        MPI_ANY_SOURCE,
                        comm_tag,
                        mp->comm,
                        &recv_req);
                TYPH_ASSERT(mpi_err == MPI_SUCCESS, ERR_MPI, TYPH_ERR_MPI,
                        "MPI_Irecv failed");
            }
        } // while (true)

        for (int i = 0; i < send_req_count; i++) {
            mpi_err = MPI_Wait(&send_req[i], &mpi_status);
            TYPH_ASSERT(mpi_err == MPI_SUCCESS, ERR_MPI, TYPH_ERR_MPI,
                    "MPI_Wait failed");
        }

        delete[] send_item_list;
        delete[] recv_item_list;

        #undef IXsil
        #undef IXril

    } // for layer_min <= layer <= layer_max

    return TYPH_SUCCESS;
}

/** @} */

} // namespace

/**
 * \param [in]  key_set_id   TODO
 * \param [out] key_set      TODO
 *
 * \returns TYPH_SUCCESS or TYPH_FAIL
 */
int
Get_Key_Set(
        int key_set_id,
        Key_Set **key_set)
{
    if ((decltype(key_sets.size())) key_set_id >= key_sets.size()) {
        return TYPH_FAIL;
    }

    if (key_set == nullptr) {
        return TYPH_FAIL;
    }

    *key_set = key_sets[key_set_id];
    return TYPH_SUCCESS;
}



/**
 * Provides a key pointer and ID matching proc and layer, given a key set and
 * instruction to search either the send or recv keys. Keys are sorted on
 * insertion by proc and then layer, and this fact is used here. This function
 * will also check for duplicate keys (which are invalid).
 *
 * \param [in]  key_set     key set to search
 * \param [in]  sendrecv    are we looking for a send or recv key?
 * \param [in]  proc        target key's proc tag
 * \param [in]  layer       target key's layer tag
 * \param [out] key_id      ID of located key
 * \param [out] key         pointer to located key
 *
 * \returns TYPH_SUCCESS if key found, TYPH_FAIL if not
 */
int
Get_Key(
        Key_Set const *key_set,
        TYPH_Sendrecv sendrecv,
        int proc,
        int layer,
        int &key_id,
        Key **key)
{
    // Initialise outputs
    key_id = -1;
    *key = nullptr;
    if (key_set == nullptr) return TYPH_FAIL;

    // Get head of key list for this key type
    Key *cur_key;
    int num_keys;
    switch (sendrecv) {
    case TYPH_SENDRECV_SEND:
        cur_key = key_set->send_keys;
        num_keys = key_set->num_send;
        break;

    case TYPH_SENDRECV_RECV:
        cur_key = key_set->recv_keys;
        num_keys = key_set->num_recv;
        break;

    default:
        TYPH_ASSERT(
                false,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Invalid sendrecv");
    }

    // Return if no key found
    if (cur_key == nullptr || num_keys <= 0) return TYPH_FAIL;

    // Find first key matching proc
    int i;
    for (i = 0; i < num_keys; i++) {
        TYPH_ASSERT(
                cur_key != nullptr,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Fewer keys than advertised");

        if (cur_key->proc == proc) {
            break;
        }

        cur_key = cur_key->next;
    }

    // Key not found. Have run off end of list
    if (i == num_keys) return TYPH_FAIL;

    // Now hunt through rest of key list for matching layer
    for (; i < num_keys; i++) {
        TYPH_ASSERT(
                cur_key != nullptr,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Fewer keys than advertised");

        if (cur_key->layer == layer) {
            break;
        }

        cur_key = cur_key->next;
    }

    // Return null if no key found
    if (cur_key == nullptr) return TYPH_FAIL;

    // Check there aren't two matching keys
    Key *tmp = cur_key->next;
    if (tmp != nullptr) {
        TYPH_ASSERT(
                (tmp->proc != cur_key->proc) ||
                    (tmp->layer != cur_key->layer),
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Duplicate key found");
    }

    key_id = i;
    *key = cur_key;
    return TYPH_SUCCESS;
}

} // namespace _TYPH_Internal

/**
 * \param [in]  key_type        key type for created keys
 * \param [in]  layer_min       min ghost layer
 * \param [in]  layer_max       max ghost layer
 * \param [in]  partition_id    ID for mesh partitioning information
 * \param [out] key_set_id      ID of the created key set
 *
 * \return TYPH_SUCCESS or TYPH_FAIL
 */
int
TYPH_Create_Key_Set(
        TYPH_Keytype key_type,
        int layer_min,
        int layer_max,
        int partition_id,
        int *key_set_id)
{
    using namespace _TYPH_Internal;

    // Get pointer to decomposition info
    Partition_Info *partition_info;
    int typh_err = Get_Partition_Info(partition_id, &partition_info);
    TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Partition info " + std::to_string(partition_id) + " not found");

    TYPH_Centring centring;
    TYPH_Auxiliary aux;
    int stride;
    int *to_glob = nullptr;
    int *to_proc = nullptr;
    int *connectivity = nullptr;

    std::unique_ptr<int[]> total_items(new int[partition_info->num_layers + 1]);
    TYPH_ALLOC_CHECK(total_items.get(), "total_items");

    switch (key_type) {
    case TYPH_KEYTYPE_CELL:
        centring = TYPH_CENTRING_CELL;
        aux      = TYPH_AUXILIARY_NONE;

        stride = 1;
        to_proc = partition_info->el_to_proc;
        to_glob = partition_info->el_loc_to_glob;

        std::copy(
                partition_info->num_el_total,
                partition_info->num_el_total + partition_info->num_layers + 1,
                total_items.get());
        break;

    case TYPH_KEYTYPE_NODE:
        centring = TYPH_CENTRING_NODE;
        aux      = TYPH_AUXILIARY_NONE;

        stride = 1;
        to_proc = partition_info->nd_to_proc;
        to_glob = partition_info->nd_loc_to_glob;

        std::copy(
                partition_info->num_nd_total,
                partition_info->num_nd_total + partition_info->num_layers + 1,
                total_items.get());
        break;

    case TYPH_KEYTYPE_CELL_CORNER:
        centring = TYPH_CENTRING_CELL;
        aux      = TYPH_AUXILIARY_NONE;

        stride = partition_info->nodes_per_el;
        to_proc = partition_info->el_to_proc;
        to_glob = partition_info->el_loc_to_glob;
        connectivity = partition_info->connectivity;

        if (connectivity == nullptr) {
            TYPH_ASSERT(
                    false,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Connectivity array must be supplied for keytype = "
                        "TYPH_CELL_CORNER");
        }

        std::copy(
                partition_info->num_el_total,
                partition_info->num_el_total + partition_info->num_layers + 1,
                total_items.get());
        break;

    default:
        TYPH_ASSERT(false, ERR_INT, TYPH_ERR_INTERNAL, "Invalid key type");
    }

    // Check limits are correct
    if (centring == TYPH_CENTRING_NODE) {
        TYPH_ASSERT(
                ((layer_min >= 0) && (layer_min <= partition_info->num_layers)),
                ERR_INT, TYPH_ERR_INTERNAL,
                "layer_min is not between 0 and num_layers");
        TYPH_ASSERT(
                ((layer_max >= 0) && (layer_max <= partition_info->num_layers)),
                ERR_INT, TYPH_ERR_INTERNAL,
                "layer_max is not between 0 and num_layers");

    } else {
        TYPH_ASSERT(
                ((layer_min >= 1) && (layer_min <= partition_info->num_layers)),
                ERR_INT, TYPH_ERR_INTERNAL,
                "layer_min is not between 1 and num_layers");
        TYPH_ASSERT(
                ((layer_max >= 1) && (layer_max <= partition_info->num_layers)),
                ERR_INT, TYPH_ERR_INTERNAL,
                "layer_max is not between 1 and num_layers");
    }

    TYPH_ASSERT((layer_min <= layer_max), ERR_INT, TYPH_ERR_INTERNAL,
            "layer_min > layer_max");

    // Create a new key set to work with
    Key_Set *key_set = nullptr;
    typh_err = New_Key_Set(
            *key_set_id,
            centring,
            aux,
            stride,
            layer_min,
            layer_max,
            partition_id,
            key_set);
    TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "New_Key_Set failed");

    // Define the key set
    typh_err = Define_Key_Set(
            key_set,
            key_type,
            layer_min,
            layer_max,
            centring,
            stride,
            to_proc,
            to_glob,
            total_items.get(),
            connectivity);
    TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Define_Key_Set failed");

    return TYPH_SUCCESS;
}
