#include <typhon.h>

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
    int typh_err;
    typh_err = TYPH_Init(nullptr);
    TYPH_FAIL_ON_ERR(typh_err)

    int nproc;
    typh_err = TYPH_Get_Size(&nproc);
    TYPH_FAIL_ON_ERR(typh_err)

    int rank;
    typh_err = TYPH_Get_Rank(&rank);
    TYPH_FAIL_ON_ERR(typh_err)

    // Check that sum reduction works by adding up ranks
    rank++;
    int res;
    typh_err = TYPH_Reduce_i(&rank, nullptr, 0, &res, TYPH_OP_SUM);
    TYPH_FAIL_ON_ERR(typh_err)

    // ... Result should be nproc*(nproc+1)/2
    int const correct = nproc * (nproc + 1) / 2;
    bool const success = res == correct;

    typh_err = TYPH_Kill(1);
    TYPH_FAIL_ON_ERR(typh_err)

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
