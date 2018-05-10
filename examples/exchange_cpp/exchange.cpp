#include <iostream>
#include <cmath>
#include <cassert>

#include <typhon.h>



// Row-major index
#define IXr(stride, i, j) ((i)*(stride)+(j))



void
check_typh_err(int typh_err, std::string func)
{
    if (typh_err != TYPH_SUCCESS) {
        std::cerr << func << " failed\n";
        std::exit(1);
    }
}



int
get_el_num(int meshw, int meshh, int x, int y)
{
    assert(x < meshw);
    assert(y < meshh);

    return y * meshw + x;
}



int
get_nd_num(int meshw, int meshh, int x, int y)
{
    assert(x <= meshw);
    assert(y <= meshh);

    return y * (meshw+1) + x;
}



int
main(int argc, const char *argv[])
{
    // Mesh dims
    int const meshw = 100;
    int const meshh = 100;


    // -------------------------------------------------------------------------
    // INITIALISATION
    // -------------------------------------------------------------------------
    int typh_err = TYPH_Init(nullptr);
    check_typh_err(typh_err, "TYPH_Init");

    int nproc, rank;

    typh_err = TYPH_Get_Size(&nproc);
    check_typh_err(typh_err, "TYPH_Get_Size");

    typh_err = TYPH_Get_Rank(&rank);
    check_typh_err(typh_err, "TYPH_Get_Rank");

    if (TYPH_Is_Master()) {
        std::cout << "nproc = " << nproc << "\n";
    }

    typh_err = TYPH_Start_Register();
    check_typh_err(typh_err, "TYPH_Start_Register");

    /*
     *  O--                --> x=meshw
     *  |  x x x x . x x x x x
     *     x x x x . x x x x x
     *     x x x x . x x x x x
     *     x x x x . x x x x x
     *     x x x x . x x x x x
     *     . . . . . . . . . .
     *     x x x x . x x x x x
     *     x x x x . x x x x x
     *     x x x x . x x x x x
     *  |  x x x x . x x x x x
     *  v  x x x x . x x x x x
     * y=meshh
     *
     * Divide mesh into equal size squares. nproc must be a power of two and
     * both mesh dims must divide evenly into sqrt(nproc)
     */
    int const sqrt_nproc = (int) std::sqrt(nproc);
    assert(sqrt_nproc * sqrt_nproc == nproc);

    assert(meshw % sqrt_nproc == 0);
    assert(meshh % sqrt_nproc == 0);

    int const meshwl = meshw / sqrt_nproc;
    int const meshhl = meshh / sqrt_nproc;

    int const xoff = (rank % sqrt_nproc) * meshwl;
    int const yoff = (rank / sqrt_nproc) * meshhl;

    // Print partitioning
    for (int ip = 0; ip < nproc; ip++) {
        if (rank == ip) {
            std::cout << "rank " << ip << "\t";
            std::cout << "(" << xoff << ", " << yoff << ", " << meshwl
                      << ", " << meshhl << ")\n";
            std::cout.flush();
        }

        typh_err = TYPH_Barrier();
        check_typh_err(typh_err, "TYPH_Barrier");
    }

    // Distribute mesh
    int const nell = meshwl * meshhl; // local # elements
    int const nelg = meshw * meshh; // global # elements
    int const nndg = (meshw+1) * (meshh+1); // global # nodes
    int const ndperel = 4; // square elements have 4 corners
    int const nlay = 1; // 1 ghost layer

    if (TYPH_Is_Master()) {
        std::cout << "nelg = " << nelg << "\n";
        std::cout << "nndg = " << nndg << "\n";
    }

    // Set local element colour (each proc just starts with the elements it will
    // ultimately control)
    int *partition = new int[nell];
    std::fill(partition, partition + nell, rank);

    // Init connectivity data
    //
    //  - global element num
    //  - element region index
    //  - element material index
    //  - element global node num 1..4
    //
    int constexpr NDAT = 7;
    int *conn_data = new int[NDAT * nell];
    int iel = 0;
    for (int i = 0; i < meshwl; i++) {
        for (int j = 0; j < meshhl; j++) {
            int const global_el_num = get_el_num(meshw, meshh, xoff+i, yoff+j);

            // number nodes counter clockwise
            //
            // 4     3
            //  x---x
            //  |   |
            //  x---x
            // 1     2
            //
            int global_nd_num[4];
            global_nd_num[0] = get_nd_num(meshw, meshh, xoff+i  , yoff+j+1);
            global_nd_num[1] = get_nd_num(meshw, meshh, xoff+i+1, yoff+j+1);
            global_nd_num[2] = get_nd_num(meshw, meshh, xoff+i+1, yoff+j  );
            global_nd_num[3] = get_nd_num(meshw, meshh, xoff+i  , yoff+j  );

            conn_data[iel * NDAT + 0] = global_el_num;
            conn_data[iel * NDAT + 1] = 0; // don't care
            conn_data[iel * NDAT + 2] = 0; // don't care
            conn_data[iel * NDAT + 3] = global_nd_num[0];
            conn_data[iel * NDAT + 4] = global_nd_num[1];
            conn_data[iel * NDAT + 5] = global_nd_num[2];
            conn_data[iel * NDAT + 6] = global_nd_num[3];
            iel++;
        }
    }

    // Call TYPH_Distribute_Mesh
    int *layer_nel = new int[nlay + 1];
    int *layer_nnd = new int[nlay + 1];
    int total_nel = 0;
    int total_nnd = 0;

    int *ellocglob = nullptr;
    int *ndlocglob = nullptr;
    int *elreg     = nullptr;
    int *elmat     = nullptr;
    int *elnd      = nullptr;
    int *elowner   = nullptr;
    int *ndowner   = nullptr;

    typh_err = TYPH_Distribute_Mesh(nell, nelg, nndg, ndperel, nlay,
            conn_data, partition, layer_nel, layer_nnd, &total_nel, &total_nnd,
            &ellocglob, &ndlocglob, &elreg, &elmat, &elnd, &elowner, &ndowner);
    check_typh_err(typh_err, "TYPH_Distribute_Mesh");

    delete[] elmat;
    delete[] elreg;
    delete[] conn_data;
    delete[] partition;

    // Init Typhon decomposition info
    // ... Calculate cumulative layer sizes
    int *layer_cum_nel = new int[nlay+1];
    int *layer_cum_nnd = new int[nlay+1];
    layer_cum_nel[0] = layer_nel[0];
    layer_cum_nnd[0] = layer_nnd[0];
    for (int ilay = 1; ilay <= nlay; ilay++) {
        layer_cum_nel[ilay] = layer_cum_nel[ilay-1] + layer_nel[ilay];
        layer_cum_nnd[ilay] = layer_cum_nnd[ilay-1] + layer_nnd[ilay];
    }

    int partition_id;
    typh_err = TYPH_Set_Partition_Info(&partition_id, TYPH_SHAPE_QUAD, nlay,
            layer_cum_nel, layer_cum_nnd, elowner, ndowner, ellocglob,
            ndlocglob, elnd, "partitioning");
    check_typh_err(typh_err, "TYPH_Set_Partition_Info");

    delete[] layer_cum_nnd;
    delete[] layer_cum_nel;

    // Create key set
    int keyset_id;
    typh_err = TYPH_Create_Key_Set(TYPH_KEYTYPE_CELL, 1, nlay, partition_id,
            &keyset_id);
    check_typh_err(typh_err, "TYPH_Create_Key_Set");

    // Add a quant
    // Set local values to our rank and ghosts to -1
    double *elquant = new double[total_nel * ndperel];
    for (int iel = 0; iel < layer_nel[0]; iel++) {
        elquant[IXr(ndperel, iel, 0)] = (double) rank;
        elquant[IXr(ndperel, iel, 1)] = (double) rank;
        elquant[IXr(ndperel, iel, 2)] = (double) rank;
        elquant[IXr(ndperel, iel, 3)] = (double) rank;
    }

    for (int iel = layer_nel[0]; iel < total_nel; iel++) {
        elquant[IXr(ndperel, iel, 0)] = -1.;
        elquant[IXr(ndperel, iel, 1)] = -1.;
        elquant[IXr(ndperel, iel, 2)] = -1.;
        elquant[IXr(ndperel, iel, 3)] = -1.;
    }

    // Register quant with Typhon
    int dims[2] = {TYPH_MESH_DIM, ndperel};
    int elquant_id;
    typh_err = TYPH_Add_Quant(&elquant_id, "elquant", TYPH_GHOSTS_ONE,
            TYPH_DATATYPE_REAL, TYPH_CENTRING_CELL, TYPH_PURE,
            TYPH_AUXILIARY_NONE, dims, 2);
    check_typh_err(typh_err, "TYPH_Add_Quant");

    dims[0] = total_nel;
    typh_err = TYPH_Set_Quant_Address(elquant_id, elquant, dims, 2);
    check_typh_err(typh_err, "TYPH_Set_Quant_Address");

    // Add a phase
    int phase_id;
    typh_err = TYPH_Add_Phase(&phase_id, "phase", TYPH_GHOSTS_ONE, TYPH_PURE,
            keyset_id, -1);
    check_typh_err(typh_err, "TYPH_Add_Phase");

    typh_err = TYPH_Add_Quant_To_Phase(phase_id, elquant_id, -1, -1, -1, -1);
    check_typh_err(typh_err, "TYPH_Add_Quant_To_Phase");

    // Initialisation complete
    typh_err = TYPH_Finish_Register();
    check_typh_err(typh_err, "TYPH_Finish_Register");


    // -------------------------------------------------------------------------
    // EXCHANGE
    // -------------------------------------------------------------------------
    // Before exchange, all ghosts are set to -1
    for (int iel = layer_nel[0]; iel < total_nel; iel++) {
        for (int j = 0; j < ndperel; j++) {
            if (elquant[IXr(ndperel, iel, j)] != -1.) {
                std::cerr << "error\n";
                std::exit(1);
            }
        }
    }

    typh_err = TYPH_Exchange(phase_id);
    check_typh_err(typh_err, "TYPH_Exchange");

    // After exchange, ghosts should contain their owner's rank
    for (int iel = layer_nel[0]; iel < total_nel; iel++) {
        for (int j = 0; j < ndperel; j++) {
            if (elquant[IXr(ndperel, iel, j)] !=
                    (double) elowner[IXr(2, iel, 0)]) {
                std::cerr << "error\n";
                std::exit(1);
            }
        }
    }

    typh_err = TYPH_Barrier();
    check_typh_err(typh_err, "TYPH_Barrier");
    if (TYPH_Is_Master()) std::cout << "exchange succeeded\n";


    // -------------------------------------------------------------------------
    // CLEANUP
    // -------------------------------------------------------------------------
    delete[] elquant;
    delete[] ndowner;
    delete[] elowner;
    delete[] elnd;
    delete[] ndlocglob;
    delete[] ellocglob;
    delete[] layer_nnd;
    delete[] layer_nel;

    typh_err = TYPH_Kill(1);
    check_typh_err(typh_err, "TYPH_Kill");

    return EXIT_SUCCESS;
}
