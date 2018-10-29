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
#ifndef TYPHON_H
#define TYPHON_H

#include <mpi.h>



#ifdef __cplusplus
extern "C" {
#endif

/**
 * \addtogroup typhon
 *
 * @{
 */

// -----------------------------------------------------------------------------
// Public Typhon API - constants
// -----------------------------------------------------------------------------
int const TYPH_MESH_DIM = -88;  //!< Used in the public API to specify quant
                                //!  dimensions

int const TYPH_SUCCESS =  0;    //!< Function returns success
int const TYPH_FAIL    = -1;    //!< Function returns failure

int const TYPH_PURE = 1;        //!< Pure or auxiliary value


// -----------------------------------------------------------------------------
// Public Typhon API - type definitions
// -----------------------------------------------------------------------------
// XXX(timrlaw): This seems rather bookleaf specific
/** Timestep information. */
typedef struct
{
    double rdt;
    int    idt;
    char   sdt[8];
    char   mdt[10];
} TYPH_Dt;

/** Supported basic data types that Typhon can communicate. */
typedef enum {
    TYPH_DATATYPE_NULL    = 0,
    TYPH_DATATYPE_REAL    = 1,
    TYPH_DATATYPE_INTEGER = 2,
    TYPH_DATATYPE_LOGICAL = 3
} TYPH_Datatype;

// XXX(timrlaw): Not yet implemented
/** TODO */
typedef enum {
    TYPH_AUXILIARY_NONE = -9099,
} TYPH_Auxiliary;

/** Supported ghost layer counts. */
typedef enum {
    TYPH_GHOSTS_ZERO  = 0,
    TYPH_GHOSTS_ONE   = 1,
    TYPH_GHOSTS_TWO   = 2,
    TYPH_GHOSTS_THREE = 3
} TYPH_Ghosts;

/** Supported cell shapes. */
typedef enum {
    TYPH_SHAPE_QUAD = 4
} TYPH_Shape;

/** Supported quantity centrings. */
typedef enum {
    TYPH_CENTRING_NODE = 2001,
    TYPH_CENTRING_CELL = 2002
} TYPH_Centring;

/** Supported key types. */
typedef enum {
    TYPH_KEYTYPE_CELL        = 1,
    TYPH_KEYTYPE_NODE        = 2,
    TYPH_KEYTYPE_CELL_CORNER = 3
} TYPH_Keytype;

/** Available reduction operations. */
typedef enum {
    TYPH_OP_SUM  = 1001,
    TYPH_OP_PROD = 1002,
    TYPH_OP_MAX  = 1003,
    TYPH_OP_MIN  = 1004,
    TYPH_OP_OR   = 1011,
    TYPH_OP_XOR  = 1012,
    TYPH_OP_AND  = 1013
} TYPH_Op;

/** Specify a send or receive operation. */
typedef enum {
    TYPH_SENDRECV_SEND = 1001,
    TYPH_SENDRECV_RECV = 1002
} TYPH_Sendrecv;


// -----------------------------------------------------------------------------
// Public Typhon API - utility functions
// -----------------------------------------------------------------------------
/** \brief Free memory allocated by Typhon calls. */
int
TYPH_Free(int *ptr);

/** \brief Convert a #TYPH_Datatype to its string representation. */
char const *
TYPH_To_String_Datatype(TYPH_Datatype datatype);

/** \brief Convert a #TYPH_Auxiliary to its string representation. */
char const *
TYPH_To_String_Auxiliary(TYPH_Auxiliary aux);

/** \brief Convert a #TYPH_Centring to its string representation. */
char const *
TYPH_To_String_Centring(TYPH_Centring centring);


// -----------------------------------------------------------------------------
// Public Typhon API - core functions
// -----------------------------------------------------------------------------
/** Initialise Typhon. */
int
TYPH_Init(MPI_Comm *comm);

int
TYPH_Init_F(MPI_Fint *comm);

/** Shutdown Typhon. */
int
TYPH_Kill(int finalise);

/** Abort Typhon in the case of an error. */
int
TYPH_Abort(int abort_code);

/** Get the current processor count. */
int
TYPH_Get_Size(int *size);

/** Get the rank of the calling processor. */
int
TYPH_Get_Rank(int *rank);

/** Return whether or not the calling processor is the master rank. */
int
TYPH_Is_Master(void);

/** Synchronise all processors. */
int
TYPH_Barrier(void);

/** Return timestamp. */
double
TYPH_Get_Time(void);

/** Set the argument to the current MPI communicator. */
int
TYPH_Set_Comm(MPI_Comm *comm);

int
TYPH_Set_Comm_F(MPI_Fint *comm);

/** Set the argument to the MPI self communicator. */
int
TYPH_Set_Comm_Self(MPI_Comm *comm);

int
TYPH_Set_Comm_Self_F(MPI_Fint *comm);


// -----------------------------------------------------------------------------
// Public Typhon API - registration functions
// -----------------------------------------------------------------------------
/** \brief Start phase and quant registration. */
int
TYPH_Start_Register(void);

/** \brief End phase and quant registration. */
int
TYPH_Finish_Register(void);

/** \brief Undergoing phase and quant registration? */
int
TYPH_Is_Registering(int *is_registering);

/** \brief Register a new communication phase. */
int
TYPH_Add_Phase(
        int *phase_id,
        char const *name,
        TYPH_Ghosts num_ghosts,
        int pure_or_aux,
        int key_set_id,  // default -1
        int ghosts_min); // default -1

/** \brief Register a new quantity. */
int
TYPH_Add_Quant(
        int *quant_id,
        char const *name,
        TYPH_Ghosts num_ghosts,
        TYPH_Datatype datatype,
        TYPH_Centring centring,
        int pure_or_aux,
        TYPH_Auxiliary aux, // default = TYPH_AUXILIARY_NONE
        int const *dims,    // default = nullptr
        int num_dims);      // default = 0


// -----------------------------------------------------------------------------
// Public Typhon API - quantity related functions
// -----------------------------------------------------------------------------
/** \brief Add a quantity to a communication phase. */
int
TYPH_Add_Quant_To_Phase(
        int phase_id,
        int quant_id,
        int recv_quant_id, // default = -1
        int key_set_id,    // default = -1
        int ghosts_min,    // default = -1
        int ghosts_max);   // default = -1

/** \brief Set the base address for the given quantity. */
int
TYPH_Set_Quant_Address(
        int quant_id,
        void *data,
        int const *dims,
        int rank);


// -----------------------------------------------------------------------------
// Public Typhon API - decomposition related functions
// -----------------------------------------------------------------------------
/** \brief Specify how the mesh is distributed amongst processors. */
int
TYPH_Set_Partition_Info(
        int *partition_id,
        TYPH_Shape el_shape,
        int num_layers,
        int const *num_el_total,
        int const *num_nd_total,
        int const *el_to_proc,
        int const *nd_to_proc,
        int const *el_loc_to_glob,
        int const *nd_loc_to_glob,
        int const *connectivity,    // default = nullptr
        char const *name);          // default = ""


// -----------------------------------------------------------------------------
// Public Typhon API - exchange related functions
// -----------------------------------------------------------------------------
/** \brief Execute the given communication phase. */
int
TYPH_Exchange(int phase_id);

/** \brief Start the given communication phase. */
int
TYPH_Start_Exchange(int phase_id);

/** \brief Finalise the given communication phase. */
int
TYPH_Finish_Exchange(int phase_id);


// -----------------------------------------------------------------------------
// Public Typhon API - key related functions
// -----------------------------------------------------------------------------
/** \brief Create a new key set. */
int
TYPH_Create_Key_Set(
        TYPH_Keytype key_type,
        int layer_min,
        int layer_max,
        int partition_id,
        int *key_set_id);

/** \brief Get the number of send/receive indices for a given key set. */
int
TYPH_Get_Key_Set_Num_Indices(
        int key_set_id,
        TYPH_Sendrecv sendrecv,
        int *num_indices);

/** \brief Get the send/receive indices for a given key set. */
int
TYPH_Get_Key_Set_Indices(
        int key_set_id,
        TYPH_Sendrecv sendrecv,
        int *indices);


// -----------------------------------------------------------------------------
// Public Typhon API - mesh distribution
// -----------------------------------------------------------------------------
/**
 * \brief Distribute the specified mesh amongst the current processors, given a
 *        partitioning.
 */
int
TYPH_Distribute_Mesh(
        int nel_local,
        int nel_glob,
        int nnd_glob,
        int nd_per_el,
        int num_layers,
        int const *conn_data,
        int const *partition,
        int *layer_nel,
        int *layer_nnd,
        int *total_nel,
        int *total_nnd,
        int **el_loc_glob,
        int **nd_loc_glob,
        int **el_region,
        int **el_material,
        int **el_nd,
        int **el_owner,
        int **nd_owner);


// -----------------------------------------------------------------------------
// Public Typhon API - collect operations
// -----------------------------------------------------------------------------
/** Perform an all-gather operation on integer data. */
int
TYPH_Gather_i(
        int *in,
        int const *dims,
        int rank,
        int *out);

/** Perform an all-gather operation on real data. */
int
TYPH_Gather_d(
        double *in,
        int const *dims,
        int rank,
        double *out);

/** Perform an all-gather operation on boolean data. */
int
TYPH_Gather_z(
        int *in,
        int const *dims,
        int rank,
        int *out);

/** Perform an all-gather operation. */
int
TYPH_Gather(
        TYPH_Datatype datatype,
        void *in,
        int const *dims,
        int rank,
        void *out);

/** Perform an all-reduce operation on integer data. */
int
TYPH_Reduce_i(
        int *in,
        int const *dims,
        int rank,
        int *out,
        TYPH_Op op);

/** Perform an all-reduce operation on real data. */
int
TYPH_Reduce_d(
        double *in,
        int const *dims,
        int rank,
        double *out,
        TYPH_Op op);

/** Perform an all-reduce operation on boolean data. */
int
TYPH_Reduce_z(
        int *in,
        int const *dims,
        int rank,
        int *out,
        TYPH_Op op);

/** Perform an all-reduce operation. */
int
TYPH_Reduce(
        TYPH_Datatype datatype,
        void *in,
        int const *dims,
        int rank,
        void *out,
        TYPH_Op op);


// -----------------------------------------------------------------------------
// Public Typhon API - serialisation related functions
// -----------------------------------------------------------------------------
/** \brief Serialise the given key set to a file. */
int
TYPH_Serialise_Key_Set(char const *filename, int key_set_id);

/** \brief Serialise the given phase to a file. */
int
TYPH_Serialise_Phase(char const *filename, int phase_id);

/** \brief Serialise all phases to a file. */
int
TYPH_Serialise_All_Phases(char const *filename);

/** \brief Serialise the given quantity to a file. */
int
TYPH_Serialise_Quant(char const *filename, int quant_id);

/** \brief Serialise all quantities to a file. */
int
TYPH_Serialise_All_Quants(char const *filename);


// -----------------------------------------------------------------------------
// Public Typhon API - timestep reduction related functions
// -----------------------------------------------------------------------------
/** \brief Initialise dt reduction. */
int
TYPH_Add_Reduce_Dt(void);

/** \brief Reduce dt amongst processors. */
int
TYPH_Reduce_Dt(TYPH_Dt *val);


/** @} */

#ifdef __cplusplus
} // extern "C"
#endif



#endif // TYPHON_H
