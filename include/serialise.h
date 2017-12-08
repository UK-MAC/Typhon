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
#ifndef TYPHON_SERIALISE_H
#define TYPHON_SERIALISE_H



namespace _TYPH_Internal {

/**
 * \addtogroup internal
 *
 * @{
 */

static constexpr int PROC_HEADER_SERIALISATION_MARKER = 0x99;
int Write_Proc_Header(std::ostream &os, int rank, int nproc);

int Serialise(std::ostream &os, Key const *key);
int Serialise(std::ostream &os, Key_Set const *key_set);
int Serialise(std::ostream &os, Phase const *phase);
int Serialise(std::ostream &os, Schedule const *schedule);
int Serialise(std::ostream &os, Quant_Info const *quant_info);
int Serialise(std::ostream &os, Schedule_Part const *schedule_part);
int Serialise(std::ostream &os, Quant const *quant);

int Get_Marker(std::istream &is, int *marker);
int Read_Proc_Header(std::istream &is, int *rank, int *nproc);

int Deserialise(std::istream &is, Key *key);
int Deserialise(std::istream &is, Key_Set *key_set);
int Deserialise(std::istream &is, Phase *phase);
int Deserialise(std::istream &is, Schedule *schedule);
int Deserialise(std::istream &is, Quant_Info *quant_info);
int Deserialise(std::istream &is, Schedule_Part *schedule_part);
int Deserialise(std::istream &is, Quant *quant);

/** @} */

} // namespace _TYPH_Internal

// -----------------------------------------------------------------------------
// Public Typhon API - serialisation related functions
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

/** \brief Serialise the given key set to a file. */
int
TYPH_Serialise_Key_Set(
        std::string filename,
        int key_set_id);

/** \brief Serialise the given phase to a file. */
int
TYPH_Serialise_Phase(
        std::string filename,
        int phase_id);

/** \brief Serialise all phases to a file. */
int
TYPH_Serialise_All_Phases(
        std::string filename);

/** \brief Serialise the given quantity to a file. */
int
TYPH_Serialise_Quant(
        std::string filename,
        int quant_id);

/** \brief Serialise all quantities to a file. */
int
TYPH_Serialise_All_Quants(
        std::string filename);

/** @} */



#endif // TYPHON_SERIALISE_H
