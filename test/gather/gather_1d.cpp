#include <typhon.h>

#include <cstdlib>
#include <numeric>
#include <vector>
#include <memory>
#include <functional>
#include <algorithm>



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

    int constexpr N = 10;
    int data[N];
    std::iota(data, data + N, 0);
    std::transform(data, data + N, data, [rank](int v) { return v * rank; });

    // Check that gather works by gathering ranks
    std::unique_ptr<int[]> res(new int[nproc * N]);
    typh_err = TYPH_Gather_i(&data[0], &N, 1, res.get());
    TYPH_FAIL_ON_ERR(typh_err)

    bool success = true;
    for (int i = 0; i < nproc; i++) {
        for (int j = 0; j < N; j++) {
            if (res[i*N + j] != j * i) {
                success = false;
                break;
            }
        }

        if (!success) break;
    }

    typh_err = TYPH_Kill(1);
    TYPH_FAIL_ON_ERR(typh_err)

    return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
