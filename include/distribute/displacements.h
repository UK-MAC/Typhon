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
#ifndef TYPHON_DISTRIBUTE_DISPLACEMENTS_H
#define TYPHON_DISTRIBUTE_DISPLACEMENTS_H

#include <numeric>

#include "typhon.h"



namespace _TYPH_Internal {

/**
 * Structure for holding counts and displacements per processor, used in MPI
 * operations.
 */
struct Displacements
{
    int *counts = nullptr;
    int *displs = nullptr;
    int nproc = 0;

    explicit Displacements(int _nproc) : nproc(_nproc)
    {
        counts = new int[nproc];
        displs = new int[nproc];

        std::fill(counts, counts + nproc, 0);
        std::fill(displs, displs + nproc, 0);
    }

    ~Displacements()
    {
        #define SAFE_DELETE(ptr) \
            if ((ptr) != nullptr) { \
                delete[] (ptr); \
                (ptr) = nullptr; \
            }

        SAFE_DELETE(counts);
        SAFE_DELETE(displs);

        #undef SAFE_DELETE
    }

    void
    zeroCounts()
    {
        std::fill(counts, counts + nproc, 0);
    }

    int
    sumCounts()
    {
        return std::accumulate(counts, counts + nproc, 0);
    }

    void
    multiplyCounts(int a)
    {
        std::transform(counts, counts + nproc, counts,
                [a](int v) { return a * v; });
    }

    void
    divideCounts(int a)
    {
        std::transform(counts, counts + nproc, counts,
                [a](int v) { return v / a; });
    }

    void
    update()
    {
        displs[0] = 0;
        for (int i = 1; i < nproc; i++) {
            displs[i] = displs[i-1] + counts[i-1];
        }
    }
};



inline int
exchangeCounts(
        Displacements &sdispls,
        Displacements &rdispls,
        MPI_Comm *comm)
{
    TYPH_ASSERT(
            sdispls.nproc == rdispls.nproc,
            ERR_INT,
            TYPH_ERR_INTERNAL,
            "mismatched nproc");

    int const mpi_err = MPI_Alltoall(
            sdispls.counts,
            1,
            MPI_INT,
            rdispls.counts,
            1,
            MPI_INT,
            *comm);

    if (mpi_err != MPI_SUCCESS) {
        TYPH_ASSERT(
                false,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Alltoall failed");
        return -1;
    }

    // recv_count
    return std::accumulate(rdispls.counts, rdispls.counts + rdispls.nproc, 0);
}



inline void
exchangeData(
        int *sdata,
        int *rdata,
        Displacements &sdispls,
        Displacements &rdispls,
        MPI_Comm *comm)
{
    int const mpi_err = MPI_Alltoallv(
            sdata,
            sdispls.counts,
            sdispls.displs,
            MPI_INT,
            rdata,
            rdispls.counts,
            rdispls.displs,
            MPI_INT,
            *comm);

    if (mpi_err != MPI_SUCCESS) {
        TYPH_ASSERT(
                false,
                ERR_MPI,
                TYPH_ERR_MPI,
                "MPI_Alltoallv failed");
    }
}

} // namespace _TYPH_Internal



#endif // TYPHON_DISTRIBUTE_DISPLACEMENTS_H
