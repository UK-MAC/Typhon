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
    typh_err = TYPH_Init();
    TYPH_FAIL_ON_ERR(typh_err)

    int nproc;
    typh_err = TYPH_Get_Size(&nproc);
    TYPH_FAIL_ON_ERR(typh_err)

    int rank;
    typh_err = TYPH_Get_Rank(&rank);
    TYPH_FAIL_ON_ERR(typh_err)

    // Check that logical and reduction works
    // ... Try with some trues and some falses
    bool input = rank % 2 == 0;
    bool res;
    typh_err = TYPH_Reduce(&input, nullptr, 0, &res, TYPH_OP_AND);
    TYPH_FAIL_ON_ERR(typh_err)
    bool success = !res;

    // ... Try with all false
    input = false;
    typh_err = TYPH_Reduce(&input, nullptr, 0, &res, TYPH_OP_AND);
    TYPH_FAIL_ON_ERR(typh_err)
    success &= !res;

    // ... Try with all true
    input = true;
    typh_err = TYPH_Reduce(&input, nullptr, 0, &res, TYPH_OP_AND);
    TYPH_FAIL_ON_ERR(typh_err)
    success &= res;

    typh_err = TYPH_Kill();
    TYPH_FAIL_ON_ERR(typh_err)

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
