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

    // Check that max reduction works
    int res;
    typh_err = TYPH_Reduce(&rank, nullptr, 0, &res, TYPH_OP_MAX);
    TYPH_FAIL_ON_ERR(typh_err)

    int const correct = nproc - 1;
    bool const success = res == correct;

    typh_err = TYPH_Kill();
    TYPH_FAIL_ON_ERR(typh_err)

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
