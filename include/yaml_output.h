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
#ifndef TYPHON_YAML_OUTPUT_H
#define TYPHON_YAML_OUTPUT_H

#ifndef TYPHON_YAMLCPP_SUPPORT
static_assert(false, "yaml-cpp support not enabled");
#endif

#include <iostream>

#include <yaml-cpp/yaml.h>



YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Key const &key);

YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Key_Set const &key_set);

YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Quant_Info const &quant_info);

YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Schedule_Part const &schedule_part);

YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Schedule const &schedule);

YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Phase const &phase);

YAML::Emitter &
operator<<(YAML::Emitter &emitter, _TYPH_Internal::Quant const &quant);



#endif // TYPHON_YAML_OUTPUT_H
