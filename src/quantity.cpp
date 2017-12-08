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

#include <algorithm>
#include <cassert>



int
TYPH_Add_Quant_To_Phase(int phase_id, int quant_id, int recv_quant_id,
        int key_set_id, int ghosts_min, int ghosts_max)
{
    using namespace _TYPH_Internal;

    // Get phase
    Phase *phase;
    int irc = TYPH_REGISTRY->Get_Phase(phase_id, &phase);

    // Get quant
    Quant *quant;
    irc |= TYPH_REGISTRY->Get_Quant(quant_id, &quant);

    if (irc != TYPH_SUCCESS) return irc;

    // Add new quant info to the end of the quant_info linked list
    Quant_Info *quant_info = nullptr;
    if (phase->quant_info == nullptr) {
        phase->quant_info = new Quant_Info();
        TYPH_ALLOC_CHECK(phase->quant_info, "phase.quant_info");
        quant_info = phase->quant_info;

    } else {
        quant_info = phase->quant_info;
        while (true) {
            if (quant_info->next == nullptr) break;
            quant_info = quant_info->next;
        }

        quant_info->next = new Quant_Info();
        TYPH_ALLOC_CHECK(quant_info->next, "quant_info.next");
        quant_info = quant_info->next;
    }

    quant_info->next = nullptr;

    quant_info->layer_min = ghosts_min != -1 ? ghosts_min : phase->layer_min;
    quant_info->layer_max = ghosts_max != -1 ? ghosts_max : phase->layer_max;

    quant_info->quant_id = quant_id;

    Quant *quant_recv = quant;
    quant_info->recv_quant_id = quant_id;

    // Get different receiving quant if required. Check they match
    if (recv_quant_id != -1 && quant_id != recv_quant_id) {
        quant_info->recv_quant_id = recv_quant_id;
        irc = TYPH_REGISTRY->Get_Quant(recv_quant_id, &quant_recv);

        TYPH_ASSERT((quant->centring == quant_recv->centring),
                ERR_USER, TYPH_ERR_INVALID_ARG, "Centrings do not match");
        TYPH_ASSERT((quant->datatype == quant_recv->datatype),
                ERR_USER, TYPH_ERR_INVALID_ARG, "Datatypes do not match");
        TYPH_ASSERT((quant->aux == quant_recv->aux),
                ERR_USER, TYPH_ERR_INVALID_ARG, "Aux IDs do not match");
        TYPH_ASSERT((quant->rank == quant_recv->rank),
                ERR_USER, TYPH_ERR_INVALID_ARG, "Ranks do not match");

        if (quant->rank > 1) {
            bool dims_match = true;
            for (int i = 0; i < quant->rank; i++) {
                if (quant->dims[i] != quant_recv->dims[i]) {
                    dims_match = false;
                    break;
                }
            }

            TYPH_ASSERT(dims_match, ERR_USER, TYPH_ERR_INVALID_ARG,
                    "Dims do not match");
        }
    }

    TYPH_ASSERT((quant_info->layer_min <= quant_recv->num_layers),
            ERR_USER, TYPH_ERR_INVALID_ARG, "layer_min > quant.num_layers");
    TYPH_ASSERT(quant_info->layer_max <= quant_recv->num_layers,
            ERR_USER, TYPH_ERR_INVALID_ARG, "layer_max > quant.num_layers");

    quant_info->key_set_id = key_set_id != -1 ? key_set_id : phase->key_set_id;

    quant_info->key_set = nullptr;
    quant_info->quant_size = 0;

    // Add Quant to Phase list
    irc = TYPH_REGISTRY->Tag_Quant(quant_id, phase_id);

    return irc;
}



int
TYPH_Set_Quant_Address(int quant_id, void *data, int const *dims, int rank)
{
    using namespace _TYPH_Internal;

    MPI_Aint addr = TYPH_NULL_ADDR;
    if (data != nullptr && dims != nullptr) {
        MPI_Get_address(data, &addr);
    }

    // Get the Quant
    Quant *quant;
    int typh_err = TYPH_REGISTRY->Get_Quant(quant_id, &quant);
    if (typh_err != TYPH_SUCCESS) return TYPH_FAIL;

    TYPH_ASSERT(
            quant != nullptr,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Quant null");

    // Set the bounds
    if (quant->upper_bound != nullptr) {
        delete[] quant->upper_bound;
        quant->upper_bound = nullptr;
    }

    if (quant->lower_bound != nullptr) {
        delete[] quant->lower_bound;
        quant->lower_bound = nullptr;
    }

    if (addr != TYPH_NULL_ADDR) {
        quant->upper_bound = new int[rank];
        quant->lower_bound = new int[rank];

        for (int i = 0; i < rank; i++) {
            quant->lower_bound[i] = 0;
            quant->upper_bound[i] = dims[i];
        }
    }

    quant->quant_address = addr;
    return TYPH_SUCCESS;
}
