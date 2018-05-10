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
    int typh_err;
    typh_err = TYPH_Init(nullptr);
    TYPH_FAIL_ON_ERR(typh_err)

    int nproc;
    typh_err = TYPH_Get_Size(&nproc);
    TYPH_FAIL_ON_ERR(typh_err)

    int rank;
    typh_err = TYPH_Get_Rank(&rank);
    TYPH_FAIL_ON_ERR(typh_err)

    // Check that logical or reduction works
    // ... Try with some trues and some falses
    int input = rank % 2 == 0 ? 1 : 0;
    int res;
    typh_err = TYPH_Reduce_z(&input, nullptr, 0, &res, TYPH_OP_OR);
    TYPH_FAIL_ON_ERR(typh_err)
    bool success = res == 1;

    // ... Try with all false
    input = 0;
    typh_err = TYPH_Reduce_z(&input, nullptr, 0, &res, TYPH_OP_OR);
    TYPH_FAIL_ON_ERR(typh_err)
    success &= res == 0;

    // ... Try with all true
    input = 1;
    typh_err = TYPH_Reduce_z(&input, nullptr, 0, &res, TYPH_OP_OR);
    TYPH_FAIL_ON_ERR(typh_err)
    success &= res == 1;

    typh_err = TYPH_Kill(1);
    TYPH_FAIL_ON_ERR(typh_err)

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
