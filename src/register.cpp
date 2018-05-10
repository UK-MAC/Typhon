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
#include <string>
#include <map>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "typhon.h"
#include "types.h"
#include "utilities.h"
#include "core.h"
#include "register.h"



namespace _TYPH_Internal {

Registry *Registry::singleton = nullptr;



void
Registry::Init_Singleton()
{
    // Check if already initialised
    if (singleton != nullptr) return;

    singleton = new Registry();
}



void
Registry::Kill_Singleton()
{
    // Check if uninitialised
    if (singleton == nullptr) return;

    delete singleton;
}



int
Registry::Get_Phase(int phase_id, Phase **phase)
{
    if ((decltype(phases.size())) phase_id >= phases.size()) {
        return TYPH_FAIL;
    }

    if (phase != nullptr) {
        *phase = &phases[phase_id];
    }

    return TYPH_SUCCESS;
}



int
Registry::Get_Quant(int quant_id, Quant **quant)
{
    if ((decltype(quants.size())) quant_id >= quants.size()) {
        return TYPH_FAIL;
    }

    if (quant != nullptr) {
        *quant = &quants[quant_id];
    }

    return TYPH_SUCCESS;
}



int
Registry::Initialise()
{
    if (initialised) return TYPH_SUCCESS;

    // Register kill routine to Typhon core
    int irc = TYPH_CORE->Kill_Insert([this]() -> int { return Kill(); });

    phases.clear();
    quants.clear();

    initialised  = true;
    return irc;
}



int
Registry::Kill()
{
    if (!initialised) return TYPH_SUCCESS;

    int irc = TYPH_SUCCESS;
    if (registration) {
        irc |= Finish_Registration();
    }

    // Deallocate all the finalized Phase & Quant array structures
    if (phases.size() > 0) {
        for (Phase &p : phases) {
            if (p.num_quants > 0) {
                delete[] p.quant_id_list;
            }
        }

        phases.clear();
    }

    quants.clear();

    registration = false;
    initialised  = false;
    return irc;
}



int
Registry::Start_Registration()
{
    if (!initialised) return TYPH_FAIL;
    if (registration) return TYPH_SUCCESS;

    registration = true;
    return TYPH_SUCCESS;
}



int
Registry::Finish_Registration()
{
    if (!initialised)  return TYPH_FAIL;
    if (!registration) return TYPH_SUCCESS;

    // Keep a tally of the number of Quants of each datatype as this will be
    // used to allocate the QuantData stores
    std::map<TYPH_Datatype, int> quants_per_type;
    quants_per_type[TYPH_DATATYPE_REAL] = 0;
    quants_per_type[TYPH_DATATYPE_INTEGER] = 0;
    quants_per_type[TYPH_DATATYPE_LOGICAL] = 0;

    for (Quant &q : quants) {
        quants_per_type[q.datatype]++;
        q.quant_data_id = quants_per_type[q.datatype] - 1;
    }

    registration = false;
    return TYPH_SUCCESS;
}



int
Registry::Add_Phase(int &phase_id, std::string name, int num_layers,
        int pure_or_aux __attribute__((unused)), int key_set_id, int layer_min)
{
    if (!registration) return TYPH_FAIL;

    phases.push_back({});
    phase_id = phases.size() - 1;
    Phase &phase = phases[phase_id];

    phase.num_layers = num_layers;
    phase.num_quants = 0;
    //!iPhase%pure = (iand(PureOrAux, TYPH_PURE) == TYPH_PURE)
    phase.is_pure = true;
    phase.name = name;

    // Set key set ID if provided
    phase.key_set_id = key_set_id >= 0 ? key_set_id : -1;

    // Set minimum ghosts if provided
    phase.layer_min = layer_min >= 0 ? layer_min : 1;

    phase.layer_max = num_layers;
    phase.is_built = false;
    phase.is_committed = false;

    return TYPH_SUCCESS;
}



int
Registry::Add_Quant(int &quant_id, std::string name, int num_layers,
        TYPH_Datatype datatype, TYPH_Centring centring,
        int pure_or_aux __attribute__((unused)), TYPH_Auxiliary aux,
        int const *dims, int rank)
{
    if (!registration) return TYPH_FAIL;

    Quant quant;

    quant.num_layers = num_layers;
    quant.datatype = datatype;
    quant.centring = centring;
    //!iQuant%pure     = (iand(PureOrAux, TYPH_PURE) == TYPH_PURE)
    quant.is_pure = true;
    quant.name = name;

    quant.mpi_datatype = MPI_Datatype_From_TYPH_Datatype(datatype);
    TYPH_ASSERT(
            quant.mpi_datatype != MPI_DATATYPE_NULL,
            ERR_USER,
            TYPH_ERR_INVALID_OP,
            "Unknown datatype");

    quant.aux = aux;

    // Set dims
    if (dims != nullptr && rank > 1) {
        quant.rank = rank;
        quant.dims = new int[quant.rank];
        TYPH_ALLOC_CHECK(quant.dims, "quant.dims");

        for (int i = 0; i < quant.rank; i++) {
            quant.dims[i] = dims[i];
            if (dims[i] == TYPH_MESH_DIM) {
                quant.mesh_dim = i;
            }
        }

    // If dims is not provided, we assume a 1D quant along the mesh
    } else {
        quant.rank = 1;
        quant.dims = new int[quant.rank];
        TYPH_ALLOC_CHECK(quant.dims, "quant.dims");
        quant.dims[0] = TYPH_MESH_DIM;
        quant.mesh_dim = 0;
    }

    quant.quant_address = TYPH_NULL_ADDR;

    quants.push_back(quant);
    quant_id = quants.size() - 1;

    return TYPH_SUCCESS;
}



int
Registry::Tag_Quant(int quant_id, int phase_id)
{
    // Ensure the quant and phase IDs are valid
    int typh_err = Get_Quant(quant_id, nullptr);
    TYPH_ASSERT_RET(
                typh_err == TYPH_SUCCESS,
                ERR_USER,
                TYPH_ERR_INVALID_ARG,
                "Invalid quant_id");

    Phase *phase;
    typh_err = Get_Phase(phase_id, &phase);
    TYPH_ASSERT_RET(
            typh_err == TYPH_SUCCESS,
            ERR_USER,
            TYPH_ERR_INVALID_ARG,
            "Invalid phase_id");

    // Check if this quant has already been tagged
    if (phase->num_quants > 0) {
        if (std::any_of(
                    phase->quant_id_list, phase->quant_id_list + phase->num_quants,
                    [quant_id]
                    (int v) {
                        return v == quant_id;
                    })) {
            return TYPH_SUCCESS;
        }
    }

    // Make a temporary copy of the quant list
    int *tmp_list = nullptr;
    if (phase->num_quants > 0) {
        tmp_list = new int[phase->num_quants];
        std::copy(phase->quant_id_list, phase->quant_id_list + phase->num_quants,
                tmp_list);
    }

    // Resize the quant list and copy the old entries back
    if (phase->quant_id_list != nullptr) delete[] phase->quant_id_list;
    phase->quant_id_list = new int[phase->num_quants + 1];

    if (phase->num_quants > 0) {
        std::copy(tmp_list, tmp_list + phase->num_quants, phase->quant_id_list);
    }

    // Add the new quant
    phase->quant_id_list[phase->num_quants] = quant_id;
    phase->num_quants++;

    return TYPH_SUCCESS;
}



Registry::Registry() : initialised(false), registration(false)
{
}



Registry::~Registry()
{
}

} // namespace _TYPH_Internal

int
TYPH_Start_Register(void)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            !(TYPH_REGISTRY != nullptr && TYPH_REGISTRY->Is_Registering()),
            ERR_USER,
            TYPH_ERR_INVALID_OP,
            "Already in registration mode");

    // Initialise the register
    if (TYPH_REGISTRY == nullptr) {
        Registry::Init_Singleton();
        int typh_err = TYPH_REGISTRY->Initialise();
        TYPH_ASSERT_RET(
                typh_err == TYPH_SUCCESS,
                ERR_INT,
                TYPH_ERR_INTERNAL,
                "Registry::Initialise failed");
    }

    // Flag that finalised Phase / Quant structures don't yet exist
    int const typh_err = TYPH_REGISTRY->Start_Registration();
    return TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Registry::Start_Registration failed");
}



int
TYPH_Finish_Register(void)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_REGISTRY != nullptr && TYPH_REGISTRY->Is_Registering()),
            ERR_USER,
            TYPH_ERR_INVALID_OP,
            "Not in registration mode");

    int const typh_err = TYPH_REGISTRY->Finish_Registration();
    return TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Registry::Finish_Registration failed");
}



int
TYPH_Is_Registering(int *is_registering)
{
    *is_registering =
        (TYPH_REGISTRY != nullptr && TYPH_REGISTRY->Is_Registering());
    return TYPH_SUCCESS;
}



int
TYPH_Add_Phase(
        int *phase_id,
        char const *cname,
        TYPH_Ghosts num_ghosts,
        int pure_or_aux __attribute__((unused)),
        int key_set_id,
        int ghosts_min)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_REGISTRY != nullptr && TYPH_REGISTRY->Is_Registering()),
            ERR_USER,
            TYPH_ERR_INVALID_OP,
            "Not in registration mode");

    std::string name;
    if (cname != nullptr) {
        name = std::string(cname);
    }

    int num_layers;
    switch (num_ghosts) {
    case TYPH_GHOSTS_ZERO:  num_layers = 0; break;
    case TYPH_GHOSTS_ONE:   num_layers = 1; break;
    case TYPH_GHOSTS_TWO:   num_layers = 2; break;
    case TYPH_GHOSTS_THREE: num_layers = 3; break;
    default:
        TYPH_ASSERT_RET(
                false,
                ERR_USER,
                TYPH_ERR_INVALID_ARG,
                "Unrecognised value for num_ghosts");
    }

    int const typh_err = TYPH_REGISTRY->Add_Phase(*phase_id, name, num_layers,
            pure_or_aux, key_set_id, ghosts_min);
    return TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Registry::Add_Phase failed");
}



int
TYPH_Add_Quant(
        int *quant_id,
        char const *cname,
        TYPH_Ghosts num_ghosts,
        TYPH_Datatype datatype,
        TYPH_Centring centring,
        int pure_or_aux __attribute__((unused)),
        TYPH_Auxiliary aux,
        int const *dims,
        int rank)
{
    using namespace _TYPH_Internal;
    TYPH_ASSERT_RET(
            (TYPH_REGISTRY != nullptr && TYPH_REGISTRY->Is_Registering()),
            ERR_USER,
            TYPH_ERR_INVALID_OP,
            "Not in registration mode");

    std::string name;
    if (cname != nullptr) {
        name = std::string(cname);
    }

    int num_layers;
    switch (num_ghosts) {
    case TYPH_GHOSTS_ZERO:  num_layers = 0; break;
    case TYPH_GHOSTS_ONE:   num_layers = 1; break;
    case TYPH_GHOSTS_TWO:   num_layers = 2; break;
    case TYPH_GHOSTS_THREE: num_layers = 3; break;
    default:
        TYPH_ASSERT_RET(
                false,
                ERR_USER,
                TYPH_ERR_INVALID_ARG,
                "Unrecognised value for num_ghosts");
    }

    int const typh_err = TYPH_REGISTRY->Add_Quant(*quant_id, name, num_layers,
            datatype, centring, pure_or_aux, aux, dims, rank);
    return TYPH_ASSERT(
            typh_err == TYPH_SUCCESS,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "Registry::Add_Quant failed");
}
