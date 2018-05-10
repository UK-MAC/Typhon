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
#include <type_traits>
#include <fstream>
#include <numeric>

#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "core.h"
#include "register.h"
#include "keys.h"
#include "serialise.h"



namespace _TYPH_Internal {

int
Write_Proc_Header(std::ostream &os, int rank, int nproc)
{
    int const marker = PROC_HEADER_SERIALISATION_MARKER;
    os.write((char const *) &marker, sizeof(marker));

    os.write((char const *) &rank,  sizeof(rank));
    os.write((char const *) &nproc, sizeof(nproc));

    return TYPH_SUCCESS;
}



int
Serialise(std::ostream &os, Key const *key)
{
    int const marker = Key::SERIALISATION_MARKER;
    os.write((char const *) &marker, sizeof(marker));

    os.write((char const *) &key->layer, sizeof(key->layer));
    os.write((char const *) &key->proc,  sizeof(key->proc));

    os.write((char const *) &key->list_len, sizeof(key->list_len));
    os.write((char const *) key->list, key->list_len * sizeof(key->list[0]));

    // TODO(timrlaw): block_lens

    return TYPH_SUCCESS;
}



int
Serialise(std::ostream &os, Key_Set const *key_set)
{
    int const marker = Key_Set::SERIALISATION_MARKER;
    os.write((char const *) &marker, sizeof(marker));

    // Meta
    os.write((char const *) &key_set->centring,  sizeof(key_set->centring));
    os.write((char const *) &key_set->aux,       sizeof(key_set->aux));
    os.write((char const *) &key_set->stride,    sizeof(key_set->stride));
    os.write((char const *) &key_set->layer_min, sizeof(key_set->layer_min));
    os.write((char const *) &key_set->layer_max, sizeof(key_set->layer_max));
    os.write((char const *) &key_set->partition_id,
            sizeof(key_set->partition_id));

    // Keys
    os.write((char const *) &key_set->num_send, sizeof(key_set->num_send));
    os.write((char const *) &key_set->num_recv, sizeof(key_set->num_recv));

    auto Serialise_Keys = [](std::ostream &os, Key const *head, int len) {
        Key const *key = head;
        for (int i = 0; i < len; i++) {
            TYPH_ASSERT(
                    key != nullptr,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Fewer keys than anticipated");

            Serialise(os, key);
            key = key->next;
        }
    };

    Serialise_Keys(os, key_set->send_keys, key_set->num_send);
    Serialise_Keys(os, key_set->recv_keys, key_set->num_recv);

    return TYPH_SUCCESS;
}



int
Serialise(std::ostream &os, Phase const *phase)
{
    int const marker = Phase::SERIALISATION_MARKER;
    os.write((char const *) &marker, sizeof(marker));

    // Name
    int const name_len = phase->name.size();
    os.write((char const *) &name_len,       sizeof(int));
    os.write((char const *) &phase->name[0], name_len * sizeof(char));

    // Pure or aux?
    unsigned char const is_pure = phase->is_pure ? 1 : 0;
    os.write((char const *) &is_pure, sizeof(is_pure));

    // Quant IDs
    os.write((char const *) &phase->num_quants, sizeof(phase->num_quants));
    if (phase->num_quants > 0) {
        os.write((char const *) phase->quant_id_list,
                phase->num_quants * sizeof(phase->quant_id_list[0]));
    }

    // Key set
    os.write((char const *) &phase->key_set_id, sizeof(phase->key_set_id));

    // Ghost layers
    os.write((char const *) &phase->num_layers, sizeof(phase->num_layers));
    os.write((char const *) &phase->layer_min, sizeof(phase->layer_min));
    os.write((char const *) &phase->layer_max, sizeof(phase->layer_max));

    // Quant info
    if (phase->num_quants > 0) {
        Quant_Info *quant_info = phase->quant_info;
        for (int i = 0; i < phase->num_quants; i++) {
            TYPH_ASSERT(
                    quant_info != nullptr,
                    ERR_INT,
                    TYPH_ERR_INTERNAL,
                    "Fewer quant info instances than expected");

            Serialise(os, quant_info);
            quant_info = quant_info->next;
        }
    }

    // Schedule
    unsigned char const is_built = phase->is_built ? 1 : 0;
    os.write((char const *) &is_built, sizeof(is_built));
    if (phase->is_built) {
        Serialise(os, phase->schedule);
    }

    // Is committed?
    unsigned char const is_committed = phase->is_committed ? 1 : 0;
    os.write((char const *) &is_committed, sizeof(is_committed));

    return TYPH_SUCCESS;
}



int
Serialise(std::ostream &os, Schedule const *schedule)
{
    int const marker = Schedule::SERIALISATION_MARKER;
    os.write((char const *) &marker, sizeof(marker));

    // Sending
    os.write((char const *) &schedule->num_send, sizeof(schedule->num_send));

    os.write((char const *) schedule->send_procs,
            schedule->num_send * sizeof(schedule->send_procs[0]));
    os.write((char const *) schedule->send_nparts,
            schedule->num_send * sizeof(schedule->send_nparts[0]));
    os.write((char const *) schedule->send_start,
            schedule->num_send * sizeof(schedule->send_start[0]));

    // Receiving
    os.write((char const *) &schedule->num_recv, sizeof(schedule->num_recv));

    os.write((char const *) schedule->recv_procs,
            schedule->num_recv * sizeof(schedule->recv_procs[0]));
    os.write((char const *) schedule->recv_nparts,
            schedule->num_recv * sizeof(schedule->recv_nparts[0]));
    os.write((char const *) schedule->recv_start,
            schedule->num_recv * sizeof(schedule->recv_start[0]));

    // Parts
    int const total_nparts =
        std::accumulate(
                schedule->send_nparts,
                schedule->send_nparts + schedule->num_send,
                0) +
        std::accumulate(
                schedule->recv_nparts,
                schedule->recv_nparts + schedule->num_recv,
                0);

    for (int i = 0; i < total_nparts; i++) {
        Serialise(os, &schedule->parts[i]);
    }

    return TYPH_SUCCESS;
}



int
Serialise(std::ostream &os, Quant_Info const *quant_info)
{
    int const marker = Quant_Info::SERIALISATION_MARKER;
    os.write((char const *) &marker, sizeof(marker));

    // IDs
    os.write((char const *) &quant_info->quant_id,
            sizeof(quant_info->quant_id));
    os.write((char const *) &quant_info->recv_quant_id,
            sizeof(quant_info->recv_quant_id));
    os.write((char const *) &quant_info->key_set_id,
            sizeof(quant_info->key_set_id));

    // Ghost layers
    os.write((char const *) &quant_info->layer_min,
            sizeof(quant_info->layer_min));
    os.write((char const *) &quant_info->layer_max,
            sizeof(quant_info->layer_max));

    // Quant meta
    os.write((char const *) &quant_info->quant_size,
            sizeof(quant_info->quant_size));
    os.write((char const *) &quant_info->nrepeat, sizeof(quant_info->nrepeat));
    os.write((char const *) &quant_info->stride,  sizeof(quant_info->stride));

    return TYPH_SUCCESS;
}



int
Serialise(std::ostream &os, Schedule_Part const *schedule_part)
{
    int const marker = Schedule_Part::SERIALISATION_MARKER;
    os.write((char const *) &marker, sizeof(marker));

    os.write((char const *) &schedule_part->quant_size,
            sizeof(schedule_part->quant_size));
    os.write((char const *) &schedule_part->nrepeat,
            sizeof(schedule_part->nrepeat));
    os.write((char const *) &schedule_part->stride,
            sizeof(schedule_part->stride));

    unsigned char const key_present = schedule_part->key != nullptr ? 1 : 0;
    os.write((char const *) &key_present, sizeof(unsigned char));
    if (schedule_part->key != nullptr) {
        Serialise(os, schedule_part->key);
    }

    return TYPH_SUCCESS;
}



int
Serialise(std::ostream &os, Quant const *quant)
{
    int const marker = Quant::SERIALISATION_MARKER;
    os.write((char const *) &marker, sizeof(marker));

    // Name
    int const name_len = quant->name.size();
    os.write((char const *) &name_len,       sizeof(int));
    os.write((char const *) &quant->name[0], name_len * sizeof(char));

    // Meta
    os.write((char const *) &quant->quant_data_id, sizeof(quant->quant_data_id));
    os.write((char const *) &quant->num_layers,    sizeof(quant->num_layers));
    os.write((char const *) &quant->centring,      sizeof(quant->centring));
    os.write((char const *) &quant->datatype,      sizeof(quant->datatype));
    os.write((char const *) &quant->aux,           sizeof(quant->aux));

    unsigned char const is_pure = quant->is_pure ? 1 : 0;
    os.write((char const *) &is_pure, sizeof(unsigned char));

    unsigned char const is_aux = quant->is_aux ? 1 : 0;
    os.write((char const *) &is_aux, sizeof(unsigned char));

    // Size
    os.write((char const *) &quant->rank,     sizeof(quant->rank));
    os.write((char const *) &quant->mesh_dim, sizeof(quant->mesh_dim));
    os.write((char const *) &quant->stride,   sizeof(quant->stride));

    if (quant->rank > 0) {
        os.write((char const *) quant->dims,
                quant->rank * sizeof(quant->dims[0]));
    }

    // TODO(timrlaw): bounds

    return TYPH_SUCCESS;
}



int
Get_Marker(std::istream &is, int *marker)
{
    // Try to read the 4-byte marker. If we hit EOF, return fail, otherwise
    // seek back to the original offset and return success
    is.read((char *) marker, sizeof(int));
    if (is.eof()) return TYPH_FAIL;

    is.seekg(-sizeof(int), is.cur);
    return TYPH_SUCCESS;
}



int
Read_Proc_Header(std::istream &is, int *rank, int *nproc)
{
    std::remove_const<decltype(PROC_HEADER_SERIALISATION_MARKER)>::type marker;
    is.read((char *) &marker, sizeof(marker));
    TYPH_ASSERT(
            marker == PROC_HEADER_SERIALISATION_MARKER,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Marker " + std::to_string(marker) + " does not match expected " +
                std::to_string(PROC_HEADER_SERIALISATION_MARKER));

    is.read((char *) rank,  sizeof(*rank));
    is.read((char *) nproc, sizeof(*nproc));

    return TYPH_SUCCESS;
}



int
Deserialise(std::istream &is, Key *key)
{
    std::remove_const<decltype(Key::SERIALISATION_MARKER)>::type marker;
    is.read((char *) &marker, sizeof(marker));
    TYPH_ASSERT(
            marker == Key::SERIALISATION_MARKER,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Marker " + std::to_string(marker) + " does not match expected " +
                std::to_string(Key::SERIALISATION_MARKER));

    is.read((char *) &key->layer, sizeof(key->layer));
    is.read((char *) &key->proc,  sizeof(key->proc));

    is.read((char *) &key->list_len,  sizeof(key->list_len));
    key->list = new int[key->list_len];
    is.read((char *) &key->list[0], key->list_len * sizeof(key->list[0]));

    // TODO(timrlaw): block_lens

    return TYPH_SUCCESS;
}



int
Deserialise(std::istream &is, Key_Set *key_set)
{
    int marker;
    is.read((char *) &marker, sizeof(marker));
    TYPH_ASSERT(
            marker == Key_Set::SERIALISATION_MARKER,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Marker " + std::to_string(marker) + " does not match expected " +
                std::to_string(Key_Set::SERIALISATION_MARKER));

    is.read((char *) &key_set->centring,  sizeof(key_set->centring));
    is.read((char *) &key_set->aux,       sizeof(key_set->aux));
    is.read((char *) &key_set->stride,    sizeof(key_set->stride));
    is.read((char *) &key_set->layer_min, sizeof(key_set->layer_min));
    is.read((char *) &key_set->layer_max, sizeof(key_set->layer_max));

    is.read((char *) &key_set->partition_id, sizeof(key_set->partition_id));

    is.read((char *) &key_set->num_send, sizeof(key_set->num_send));
    is.read((char *) &key_set->num_recv, sizeof(key_set->num_recv));

    auto Deserialise_Keys = [](std::istream &is, int len) -> Key * {
        Key *head = nullptr;
        if (len > 0) {
            Key *key = new Key();
            head = key;
            for (int i = 0; i < len; i++) {
                Deserialise(is, key);

                if (i + 1 < len) {
                    Key *next = new Key();
                    key->next = next;
                    next->prev = key;

                    key = next;
                }
            }
        }

        return head;
    };

    key_set->send_keys = Deserialise_Keys(is, key_set->num_send);
    key_set->recv_keys = Deserialise_Keys(is, key_set->num_recv);

    return TYPH_SUCCESS;
}



int
Deserialise(std::istream &is, Phase *phase)
{
    int marker;
    is.read((char *) &marker, sizeof(marker));
    TYPH_ASSERT(
            marker == Phase::SERIALISATION_MARKER,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Marker " + std::to_string(marker) + " does not match expected " +
                std::to_string(Phase::SERIALISATION_MARKER));

    // Name
    int name_len;
    is.read((char *) &name_len, sizeof(int));

    char *name = new char[name_len + 1];
    is.read((char *) name, name_len * sizeof(char));
    name[name_len] = '\0';
    phase->name = std::string(name);
    delete[] name;

    // Pure or aux?
    unsigned char is_pure;
    is.read((char *) &is_pure, sizeof(unsigned char));
    phase->is_pure = (decltype(phase->is_pure)) is_pure;

    // Quant IDs
    is.read((char *) &phase->num_quants, sizeof(phase->num_quants));
    if (phase->num_quants > 0) {
        phase->quant_id_list = new int[phase->num_quants];
        is.read((char *) phase->quant_id_list,
                phase->num_quants * sizeof(phase->quant_id_list[0]));
    } else {
        phase->quant_id_list = nullptr;
    }

    // Key set
    is.read((char *) &phase->key_set_id, sizeof(phase->key_set_id));

    // Ghost layers
    is.read((char *) &phase->num_layers, sizeof(phase->num_layers));
    is.read((char *) &phase->layer_min, sizeof(phase->layer_min));
    is.read((char *) &phase->layer_max, sizeof(phase->num_layers));

    // Quant info
    if (phase->num_quants > 0) {
        Quant_Info *quant_info = new Quant_Info();
        phase->quant_info = quant_info;
        for (int i = 0; i < phase->num_quants; i++) {
            Deserialise(is, quant_info);

            if (i + 1 < phase->num_quants) {
                Quant_Info *next = new Quant_Info();
                quant_info->next = next;

                quant_info = next;
            }
        }

    } else {
        phase->quant_info = nullptr;
    }

    // Schedule
    unsigned char is_built;
    is.read((char *) &is_built, sizeof(unsigned char));
    phase->is_built = (decltype(phase->is_built)) is_built;
    if (phase->is_built) {
        phase->schedule = new Schedule();
        Deserialise(is, phase->schedule);
    } else {
        phase->schedule = nullptr;
    }

    // Flags
    unsigned char is_committed;
    is.read((char *) &is_committed, sizeof(unsigned char));
    phase->is_committed = (decltype(phase->is_committed)) is_committed;

    return TYPH_SUCCESS;
}



int
Deserialise(std::istream &is, Schedule *schedule)
{
    int marker;
    is.read((char *) &marker, sizeof(marker));
    TYPH_ASSERT(
            marker == Schedule::SERIALISATION_MARKER,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Marker " + std::to_string(marker) + " does not match expected " +
                std::to_string(Schedule::SERIALISATION_MARKER));

    // Sending
    is.read((char *) &schedule->num_send, sizeof(schedule->num_send));

    if (schedule->num_send > 0) {
        schedule->send_procs  = new int[schedule->num_send];
        schedule->send_nparts = new int[schedule->num_send];
        schedule->send_start  = new int[schedule->num_send];

        is.read((char *) schedule->send_procs,
                schedule->num_send * sizeof(schedule->send_procs[0]));
        is.read((char *) schedule->send_nparts,
                schedule->num_send * sizeof(schedule->send_nparts[0]));
        is.read((char *) schedule->send_start,
                schedule->num_send * sizeof(schedule->send_start[0]));
    }

    // Receiving
    is.read((char *) &schedule->num_recv, sizeof(schedule->num_recv));

    if (schedule->num_recv > 0) {
        schedule->recv_procs  = new int[schedule->num_recv];
        schedule->recv_nparts = new int[schedule->num_recv];
        schedule->recv_start  = new int[schedule->num_recv];

        is.read((char *) schedule->recv_procs,
                schedule->num_recv * sizeof(schedule->recv_procs[0]));
        is.read((char *) schedule->recv_nparts,
                schedule->num_recv * sizeof(schedule->recv_nparts[0]));
        is.read((char *) schedule->recv_start,
                schedule->num_recv * sizeof(schedule->recv_start[0]));
    }

    // Parts
    int const total_nparts =
        std::accumulate(
                schedule->send_nparts,
                schedule->send_nparts + schedule->num_send,
                0) +
        std::accumulate(
                schedule->recv_nparts,
                schedule->recv_nparts + schedule->num_recv,
                0);

    if (total_nparts > 0) {
        schedule->parts = new Schedule_Part[total_nparts];
        for (int i = 0; i < total_nparts; i++) {
            Deserialise(is, &schedule->parts[i]);
        }
    }

    return TYPH_SUCCESS;
}



int
Deserialise(std::istream &is, Quant_Info *quant_info)
{
    int marker;
    is.read((char *) &marker, sizeof(marker));
    TYPH_ASSERT(
            marker == Quant_Info::SERIALISATION_MARKER,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Marker " + std::to_string(marker) + " does not match expected " +
                std::to_string(Quant_Info::SERIALISATION_MARKER));

    // IDs
    is.read((char *) &quant_info->quant_id, sizeof(quant_info->quant_id));
    is.read((char *) &quant_info->recv_quant_id,
            sizeof(quant_info->recv_quant_id));
    is.read((char *) &quant_info->key_set_id, sizeof(quant_info->key_set_id));

    // Ghost layers
    is.read((char *) &quant_info->layer_min, sizeof(quant_info->layer_min));
    is.read((char *) &quant_info->layer_max, sizeof(quant_info->layer_max));

    // Quant meta
    is.read((char *) &quant_info->quant_size, sizeof(quant_info->quant_size));
    is.read((char *) &quant_info->nrepeat,    sizeof(quant_info->nrepeat));
    is.read((char *) &quant_info->stride,     sizeof(quant_info->stride));

    return TYPH_SUCCESS;
}



int
Deserialise(std::istream &is, Schedule_Part *schedule_part)
{
    int marker;
    is.read((char *) &marker, sizeof(marker));
    TYPH_ASSERT(
            marker == Schedule_Part::SERIALISATION_MARKER,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Marker " + std::to_string(marker) + " does not match expected " +
                std::to_string(Schedule_Part::SERIALISATION_MARKER));

    is.read((char *) &schedule_part->quant_size,
            sizeof(schedule_part->quant_size));
    is.read((char *) &schedule_part->nrepeat, sizeof(schedule_part->nrepeat));
    is.read((char *) &schedule_part->stride, sizeof(schedule_part->stride));

    unsigned char key_present;
    is.read((char *) &key_present, sizeof(unsigned char));
    if (key_present) {
        schedule_part->key = new Key();
        Deserialise(is, schedule_part->key);
    }

    return TYPH_SUCCESS;
}



int
Deserialise(std::istream &is, Quant *quant)
{
    int marker;
    is.read((char *) &marker, sizeof(marker));
    TYPH_ASSERT(
            marker == Quant::SERIALISATION_MARKER,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Marker " + std::to_string(marker) + " does not match expected " +
                std::to_string(Quant::SERIALISATION_MARKER));

    // Name
    int name_len;
    is.read((char *) &name_len, sizeof(int));

    char *name = new char[name_len + 1];
    is.read((char *) name, name_len * sizeof(char));
    name[name_len] = '\0';
    quant->name = std::string(name);
    delete[] name;

    // Meta
    is.read((char *) &quant->quant_data_id, sizeof(quant->quant_data_id));
    is.read((char *) &quant->num_layers,    sizeof(quant->num_layers));
    is.read((char *) &quant->centring,      sizeof(quant->centring));
    is.read((char *) &quant->datatype,      sizeof(quant->datatype));
    is.read((char *) &quant->aux,           sizeof(quant->aux));

    unsigned char is_pure;
    is.read((char *) &is_pure, sizeof(unsigned char));
    quant->is_pure = (bool) is_pure;

    unsigned char is_aux;
    is.read((char *) &is_aux, sizeof(unsigned char));
    quant->is_aux = (bool) is_aux;

    // Size
    is.read((char *) &quant->rank,     sizeof(quant->rank));
    is.read((char *) &quant->mesh_dim, sizeof(quant->mesh_dim));
    is.read((char *) &quant->stride,   sizeof(quant->stride));

    if (quant->rank > 0) {
        quant->dims = new int[quant->rank];
        is.read((char *) quant->dims, quant->rank * sizeof(quant->dims[0]));
    }

    // TODO(timrlaw): bounds

    return TYPH_SUCCESS;
}

} // namespace _TYPH_Internal

int
TYPH_Serialise_Key_Set(char const *cfilename, int key_set_id)
{
    using namespace _TYPH_Internal;

    std::string filename;
    if (cfilename != nullptr) {
        filename = std::string(cfilename);
    }

    Key_Set *key_set;
    int typh_err = Get_Key_Set(key_set_id, &key_set);
    TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Key set " + std::to_string(key_set_id) + " not found");

    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    typh_err = TYPH_SUCCESS;
    for (int ip = 0; ip < mp->size; ip++) {
        if (ip == mp->rank) {
            auto flags = std::ios::binary;
            if (ip != 0) flags |= std::ios::app;

            std::ofstream of(filename.c_str(), flags);

            typh_err |= Write_Proc_Header(of, mp->rank, mp->size);
            typh_err |= Serialise(of, key_set);

            of.close();
        }
        TYPH_Barrier();
    }

    return typh_err;
}



int
TYPH_Serialise_Phase(char const *cfilename, int phase_id)
{
    using namespace _TYPH_Internal;

    std::string filename;
    if (cfilename != nullptr) {
        filename = std::string(cfilename);
    }

    Phase *phase;
    int typh_err = TYPH_REGISTRY->Get_Phase(phase_id, &phase);
    TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Phase " + std::to_string(phase_id) + " not found");

    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    typh_err = TYPH_SUCCESS;
    for (int ip = 0; ip < mp->size; ip++) {
        if (ip == mp->rank) {
            auto flags = std::ios::binary;
            if (ip != 0) flags |= std::ios::app;

            std::ofstream of(filename.c_str(), flags);

            typh_err |= Write_Proc_Header(of, mp->rank, mp->size);
            typh_err |= Serialise(of, phase);

            of.close();
        }
        TYPH_Barrier();
    }

    return TYPH_SUCCESS;
}



int
TYPH_Serialise_All_Phases(char const *cfilename)
{
    using namespace _TYPH_Internal;

    std::string filename;
    if (cfilename != nullptr) {
        filename = std::string(cfilename);
    }

    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    int typh_err = TYPH_SUCCESS;
    for (int ip = 0; ip < mp->size; ip++) {
        if (ip == mp->rank) {
            auto flags = std::ios::binary;
            if (ip != 0) flags |= std::ios::app;

            std::ofstream of(filename.c_str(), flags);

            typh_err |= Write_Proc_Header(of, mp->rank, mp->size);

            for (int i = 0; i < TYPH_REGISTRY->Get_Num_Phases(); i++) {
                Phase *phase;
                int typh_err = TYPH_REGISTRY->Get_Phase(i, &phase);
                TYPH_ASSERT(
                        typh_err == TYPH_SUCCESS,
                        ERR_USER,
                        TYPH_ERR_INVALID_ARG,
                        "Phase " + std::to_string(i) + " not found");

                typh_err |= Serialise(of, phase);
            }

            of.close();
        }
        TYPH_Barrier();
    }

    return TYPH_SUCCESS;
}



int
TYPH_Serialise_Quant(char const *cfilename, int quant_id)
{
    using namespace _TYPH_Internal;

    std::string filename;
    if (cfilename != nullptr) {
        filename = std::string(cfilename);
    }

    Quant *quant;
    int typh_err = TYPH_REGISTRY->Get_Quant(quant_id, &quant);
    TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Quant " + std::to_string(quant_id) + " not found");

    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    typh_err = TYPH_SUCCESS;
    for (int ip = 0; ip < mp->size; ip++) {
        if (ip == mp->rank) {
            auto flags = std::ios::binary;
            if (ip != 0) flags |= std::ios::app;

            std::ofstream of(filename.c_str(), flags);

            typh_err |= Write_Proc_Header(of, mp->rank, mp->size);
            typh_err |= Serialise(of, quant);

            of.close();
        }
        TYPH_Barrier();
    }

    return TYPH_SUCCESS;
}



int
TYPH_Serialise_All_Quants(char const *cfilename)
{
    using namespace _TYPH_Internal;

    std::string filename;
    if (cfilename != nullptr) {
        filename = std::string(cfilename);
    }

    MPI_Runtime const *mp = TYPH_CORE->Get_MPI_Runtime();

    int typh_err = TYPH_SUCCESS;
    for (int ip = 0; ip < mp->size; ip++) {
        if (ip == mp->rank) {
            auto flags = std::ios::binary;
            if (ip != 0) flags |= std::ios::app;

            std::ofstream of(filename.c_str(), flags);

            typh_err |= Write_Proc_Header(of, mp->rank, mp->size);

            for (int i = 0; i < TYPH_REGISTRY->Get_Num_Quants(); i++) {
                Quant *quant;
                int typh_err = TYPH_REGISTRY->Get_Quant(i, &quant);
                TYPH_ASSERT(
                        typh_err == TYPH_SUCCESS,
                        ERR_USER,
                        TYPH_ERR_INVALID_ARG,
                        "Quant " + std::to_string(i) + " not found");

                typh_err |= Serialise(of, quant);
            }

            of.close();
        }
        TYPH_Barrier();
    }

    return TYPH_SUCCESS;
}
