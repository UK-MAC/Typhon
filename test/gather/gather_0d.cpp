#include <typhon.h>

#include <numeric>
#include <vector>
#include <memory>



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

    // Check that gather works by gathering ranks
    std::unique_ptr<int[]> ranks(new int[nproc]);
    typh_err = TYPH_Gather(&rank, nullptr, 0, ranks.get());
    TYPH_FAIL_ON_ERR(typh_err)

    // ... Ranks should be ordered
    std::vector<int> correct;
    std::iota(correct.begin(), correct.end(), 0);
    bool const success = std::equal(correct.begin(), correct.end(), ranks.get());

    typh_err = TYPH_Kill();
    TYPH_FAIL_ON_ERR(typh_err)

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
