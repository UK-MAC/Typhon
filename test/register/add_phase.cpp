#include "typhon.h"
#include "types.h"
#include "register.h"

#include <cstdlib>



#define TYPH_FAIL_ON_ERR(terr) \
    if ((terr) != TYPH_SUCCESS) { \
        TYPH_Abort(-1); \
    }



int
main(
        int argc           __attribute__((unused)),
        char const *argv[] __attribute__((unused)))
{
    using namespace _TYPH_Internal;

    int typh_err;
    typh_err = TYPH_Init(nullptr);
    TYPH_FAIL_ON_ERR(typh_err)

    typh_err = TYPH_Start_Register();
    TYPH_FAIL_ON_ERR(typh_err)

    // Test adding a phase
    std::string            name = "test_phase";
    TYPH_Ghosts      num_ghosts = TYPH_GHOSTS_TWO;
    int             pure_or_aux = TYPH_PURE;

    int phase_id;
    typh_err = TYPH_Add_Phase(&phase_id, name.c_str(), num_ghosts, pure_or_aux,
            -1, -1);
    TYPH_FAIL_ON_ERR(typh_err)

    typh_err = TYPH_Finish_Register();
    TYPH_FAIL_ON_ERR(typh_err)

    // Should be able to retrieve the phase from the registry
    _TYPH_Internal::Phase *phase = nullptr;
    typh_err = TYPH_REGISTRY->Get_Phase(phase_id, &phase);

    bool success = false;
    if (typh_err == TYPH_SUCCESS && phase != nullptr) {
        success = true;

        // All the phases's fields should match up with what we input.
        success &= phase->name       == name;
        success &= phase->num_layers == 2;
        success &= phase->is_pure    == true;
    }

    typh_err = TYPH_Kill(1);
    TYPH_FAIL_ON_ERR(typh_err)

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
