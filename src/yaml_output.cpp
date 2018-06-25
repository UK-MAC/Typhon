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
#include "yaml_output.h"
#include "utilities.h"
#include "serialise.h"

#include <numeric>



YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Key const &key)
{
    using namespace YAML;
    emitter << BeginMap;
    emitter << Key << "type" << Value << "key";

    emitter << Key << "layer" << Value << key.layer;
    emitter << Key << "proc"  << Value << key.proc;

    emitter << Key << "list" << Value << Flow << BeginSeq;
    for (int i = 0; i < key.list_len; i++) {
        emitter << key.list[i];
    }
    emitter << EndSeq;

    emitter << EndMap;
    return emitter;
}



YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Key_Set const &key_set)
{
    using namespace YAML;
    emitter << BeginMap;
    emitter << Key << "type" << Value << "keyset";

    emitter << Key << "centring"     << Value
            << TYPH_To_String_Centring(key_set.centring);
    emitter << Key << "aux"          << Value
            << TYPH_To_String_Auxiliary(key_set.aux);
    emitter << Key << "stride"       << Value << key_set.stride;
    emitter << Key << "layer_min"    << Value << key_set.layer_min;
    emitter << Key << "layer_max"    << Value << key_set.layer_max;
    emitter << Key << "partition_id" << Value << key_set.partition_id;
    emitter << Key << "num_send"     << Value << key_set.num_send;
    emitter << Key << "num_recv"     << Value << key_set.num_recv;

    emitter << Key << "send_keys";
    emitter << Value << BeginSeq;
    for (int i = 0; i < key_set.num_send; i++) {
        emitter << key_set.send_keys[i];
    }
    emitter << EndSeq;

    emitter << Key << "recv_keys";
    emitter << Value << BeginSeq;
    for (int i = 0; i < key_set.num_recv; i++) {
        emitter << key_set.recv_keys[i];
    }
    emitter << EndSeq;

    emitter << EndMap;
    return emitter;
}



YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Quant_Info const &quant_info)
{
    using namespace YAML;
    emitter << BeginMap;
    emitter << Key << "type" << Value << "quantinfo";

    emitter << Key << "quant_id"      << Value << quant_info.quant_id;
    emitter << Key << "recv_quant_id" << Value << quant_info.recv_quant_id;
    emitter << Key << "key_set_id"    << Value << quant_info.key_set_id;
    emitter << Key << "layer_min"     << Value << quant_info.layer_min;
    emitter << Key << "layer_max"     << Value << quant_info.layer_max;
    emitter << Key << "quant_size"    << Value << quant_info.quant_size;
    emitter << Key << "nrepeat"       << Value << quant_info.nrepeat;
    emitter << Key << "stride"        << Value << quant_info.stride;

    emitter << EndMap;
    return emitter;
}



YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Schedule_Part const &schedule_part)
{
    using namespace YAML;
    emitter << BeginMap;
    emitter << Key << "type" << Value << "schedulepart";

    emitter << Key << "quant_size" << Value << schedule_part.quant_size;
    emitter << Key << "nrepeat"    << Value << schedule_part.nrepeat;
    emitter << Key << "stride"     << Value << schedule_part.stride;
    emitter << Key << "key"        << Value << *schedule_part.key;

    emitter << EndMap;
    return emitter;
}



YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Schedule const &schedule)
{
    using namespace YAML;
    emitter << BeginMap;
    emitter << Key << "type" << Value << "schedule";

    emitter << Key << "num_send" << Value << schedule.num_send;
    emitter << Key << "num_recv" << Value << schedule.num_recv;

    emitter << Key << "send_procs" << Value << Flow << BeginSeq;
    for (int i = 0; i < schedule.num_send; i++) {
        emitter << schedule.send_procs[i];
    }
    emitter << EndSeq;

    emitter << Key << "send_nparts" << Value << Flow << BeginSeq;
    for (int i = 0; i < schedule.num_send; i++) {
        emitter << schedule.send_nparts[i];
    }
    emitter << EndSeq;

    emitter << Key << "send_start" << Value << Flow << BeginSeq;
    for (int i = 0; i < schedule.num_send; i++) {
        emitter << schedule.send_start[i];
    }
    emitter << EndSeq;

    emitter << Key << "recv_procs" << Value << Flow << BeginSeq;
    for (int i = 0; i < schedule.num_recv; i++) {
        emitter << schedule.recv_procs[i];
    }
    emitter << EndSeq;

    emitter << Key << "recv_nparts" << Value << Flow << BeginSeq;
    for (int i = 0; i < schedule.num_recv; i++) {
        emitter << schedule.recv_nparts[i];
    }
    emitter << EndSeq;

    emitter << Key << "recv_start" << Value << Flow << BeginSeq;
    for (int i = 0; i < schedule.num_recv; i++) {
        emitter << schedule.recv_start[i];
    }
    emitter << EndSeq;

    int const total_nparts =
        std::accumulate(
                schedule.send_nparts,
                schedule.send_nparts + schedule.num_send,
                0) +
        std::accumulate(
                schedule.recv_nparts,
                schedule.recv_nparts + schedule.num_recv,
                0);

    emitter << Key << "parts" << Value << BeginSeq;
    for (int i = 0; i < total_nparts; i++) {
        emitter << schedule.parts[i];
    }
    emitter << EndSeq;

    emitter << EndMap;
    return emitter;
}



YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Phase const &phase)
{
    using namespace YAML;
    emitter << BeginMap;
    emitter << Key << "type" << Value << "phase";

    emitter << Key << "name"         << Value << phase.name;
    emitter << Key << "is_pure"      << Value << phase.is_pure;
    emitter << Key << "key_set_id"   << Value << phase.key_set_id;
    emitter << Key << "num_layers"   << Value << phase.num_layers;
    emitter << Key << "layer_min"    << Value << phase.layer_min;
    emitter << Key << "layer_max"    << Value << phase.layer_max;
    emitter << Key << "is_built"     << Value << phase.is_built;
    emitter << Key << "is_committed" << Value << phase.is_committed;

    emitter << Key << "quant_id_list" << Value << Flow << BeginSeq;
    for (int i = 0; i < phase.num_quants; i++) emitter << phase.quant_id_list[i];
    emitter << EndSeq;

    emitter << Key << "quant_info" << Value << BeginSeq;
    auto const *quant_info = phase.quant_info;
    for (int i = 0; i < phase.num_quants; i++) {
        emitter << *quant_info;
        quant_info = quant_info->next;
    }
    emitter << EndSeq;

    emitter << Key << "schedule" << Value;
    if (phase.is_built) emitter << *phase.schedule;
    else                emitter << "~";

    emitter << EndMap;
    return emitter;
}



YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Quant const &quant)
{
    using namespace YAML;
    emitter << BeginMap;
    emitter << Key << "type" << Value << "quant";

    emitter << Key << "name"          << Value << quant.name;
    emitter << Key << "quant_data_id" << Value << quant.quant_data_id;
    emitter << Key << "centring"      << Value
            << TYPH_To_String_Centring(quant.centring);
    emitter << Key << "datatype"      << Value
            << TYPH_To_String_Datatype(quant.datatype);
    emitter << Key << "aux"           << Value
            << TYPH_To_String_Auxiliary(quant.aux);
    emitter << Key << "is_pure"       << Value << quant.is_pure;
    emitter << Key << "is_aux"        << Value << quant.is_aux;
    emitter << Key << "rank"          << Value << quant.rank;
    emitter << Key << "mesh_dim"      << Value << quant.mesh_dim;
    emitter << Key << "stride"        << Value << quant.stride;

    emitter << Key << "dims" << Value << Flow << BeginSeq;
    for (int i = 0; i < quant.rank; i++) {
        if (quant.dims[i] == TYPH_MESH_DIM) {
            emitter << "TYPH_MESH_DIM";
        } else {
            emitter << quant.dims[i];
        }
    }
    emitter << EndSeq;

    emitter << EndMap;
    return emitter;
}

//} // namespace _TYPH_Internal
