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
#ifndef TYPHON_REGISTER_H
#define TYPHON_REGISTER_H

#include <vector>

#include "typhon.h"



namespace _TYPH_Internal {

/**
 * \addtogroup internal
 *
 * @{
 */

/** Singleton keeping track of quants and phases. */
class Registry
{
public:
    static void Init_Singleton();
    static Registry *Get_Singleton() { return singleton; }
    static void Kill_Singleton();

    bool Is_Initialised() const { return initialised; }
    bool Is_Registering() const { return registration; }

    int Get_Phase(int phase_id, Phase **phase);
    int Get_Quant(int quant_id, Quant **quant);

    int Get_Num_Phases() const { return phases.size(); }
    int Get_Num_Quants() const { return quants.size(); }

    int Initialise();
    int Kill();

    int Start_Registration();
    int Finish_Registration();

    int Add_Phase(int &phase_id, std::string name, int num_layers,
            int pure_or_aux, int key_set_id, int layer_min);
    int Add_Quant(int &quant_id, std::string name, int num_layers,
            TYPH_Datatype datatype, TYPH_Centring centring, int pure_or_aux,
            TYPH_Auxiliary aux, int const *dims, int rank);

    int Tag_Quant(int quant_id, int phase_id);

private:
    static Registry *singleton;

    explicit Registry();    // Construction only allowed through Get_Singleton()
    ~Registry();            // Destruction only allowed through Kill_Singleton()

    Registry(Registry const &rhs);              // unimplemented copy-constructor
    Registry &operator=(Registry const &rhs);   // unimplemented operator=

    bool initialised;       // Is the registry initialised
    bool registration;      // Are we in registration mode

    std::vector<Phase> phases;
    std::vector<Quant> quants;
};

#define TYPH_REGISTRY _TYPH_Internal::Registry::Get_Singleton()

/** @} */

} // namespace _TYPH_Internal



#endif // TYPHON_REGISTER_H
