#include <typhon.h>



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
    typh_err = TYPH_Init();
    TYPH_FAIL_ON_ERR(typh_err)

    typh_err = TYPH_Start_Register();
    TYPH_FAIL_ON_ERR(typh_err)

    // Test adding a quant
    std::string            name = "test_quant";
    TYPH_Ghosts      num_ghosts = TYPH_GHOSTS_TWO;
    TYPH_Datatype      datatype = TYPH_DATATYPE_REAL;
    TYPH_Centring      centring = TYPH_CENTRING_NODE;
    int             pure_or_aux = TYPH_PURE;
    TYPH_Auxiliary          aux = TYPH_AUXILIARY_NONE;

    int const N = 10;
    int const dims[2] = {TYPH_MESH_DIM, N};
    int const rank = 2;

    int quant_id;
    typh_err = TYPH_Add_Quant(quant_id, name, num_ghosts, datatype, centring,
            pure_or_aux, aux, dims, rank);
    TYPH_FAIL_ON_ERR(typh_err)

    typh_err = TYPH_Finish_Register();
    TYPH_FAIL_ON_ERR(typh_err)

    // Should be able to retrieve the quant from the registry
    _TYPH_Internal::Quant *quant = nullptr;
    typh_err = TYPH_REGISTRY->Get_Quant(quant_id, &quant);

    bool success = false;
    if (typh_err == TYPH_SUCCESS && quant != nullptr) {
        success = true;

        // All the quant's fields should match up with what we input.
        success &= quant->name       == name;
        success &= quant->num_layers == 2;
        success &= quant->centring   == centring;
        success &= quant->datatype   == datatype;
        success &= quant->aux        == aux;
        success &= quant->is_pure    == true;
        success &= quant->is_aux     == false;
        success &= quant->rank       == rank;
        success &= quant->mesh_dim   == 0;
        success &= quant->dims[0]    == TYPH_MESH_DIM;
        success &= quant->dims[1]    == N;
    }

    typh_err = TYPH_Kill();
    TYPH_FAIL_ON_ERR(typh_err)

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
