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
#ifndef TYPHON_TYPES_H
#define TYPHON_TYPES_H

#include <string>
#include <mpi.h>
#include <iostream>
#include <algorithm>
#include <cassert>



// -----------------------------------------------------------------------------
// Public API type definitions
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

// XXX(timrlaw): This seems rather bookleaf specific
/** Timestep information. */
struct TYPH_Dt
{
    double rdt;
    int    idt;
    char   sdt[8];
    char   mdt[10];
};



/** Supported basic data types that Typhon can communicate. */
enum TYPH_Datatype : int {
    TYPH_DATATYPE_NULL    = 0,
    TYPH_DATATYPE_REAL    = 1,
    TYPH_DATATYPE_INTEGER = 2,
    TYPH_DATATYPE_LOGICAL = 3
};



// XXX(timrlaw): Not yet implemented
/** TODO */
enum TYPH_Auxiliary : int {
    TYPH_AUXILIARY_NONE = -9099,
};



/** Supported ghost layer counts. */
enum TYPH_Ghosts : int {
    TYPH_GHOSTS_ZERO,
    TYPH_GHOSTS_ONE,
    TYPH_GHOSTS_TWO,
    TYPH_GHOSTS_THREE

    // XXX(timrlaw): What are these for?
    //TYPH_GHOSTS_NULL = TYPH_NULL,
    //TYPH_GHOSTS_MAX  = 4
};



/** Supported cell shapes. */
enum TYPH_Shape : int {
    TYPH_SHAPE_QUAD = 4
};



/** Supported quantity centrings. */
enum TYPH_Centring : int {
    TYPH_CENTRING_NODE = 2001,
    TYPH_CENTRING_CELL = 2002
};



/** Supported key types. */
enum TYPH_Keytype : int {
    TYPH_KEYTYPE_CELL        = 1,
    TYPH_KEYTYPE_NODE        = 2,
    TYPH_KEYTYPE_CELL_CORNER = 3
};

/** @} */



// -----------------------------------------------------------------------------
// Public API constants
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

int constexpr TYPH_MESH_DIM = -88;  //!< Used in the public API to specify quant
                                    //!  dimensions

int constexpr TYPH_SUCCESS =  0;    //!< Function returns success
int constexpr TYPH_FAIL    = -1;    //!< Function returns failure

int constexpr TYPH_PURE = 1;        //!< Pure or auxiliary value

/** @} */



namespace _TYPH_Internal {

// -----------------------------------------------------------------------------
// Internal type definitions
// -----------------------------------------------------------------------------
/**
 * \addtogroup internal
 *
 * @{
 */

/** Error codes. */
enum TYPH_Error : int {
    TYPH_ERR_USER          = -101,
    TYPH_ERR_MEMORY        = -102,
    TYPH_ERR_MPI           = -103,
    TYPH_ERR_INTERNAL      = -104,
    TYPH_ERR_APPLICATION   = -105,
    TYPH_ERR_UNINITIALISED = -110,
    TYPH_ERR_INVALID_ARG   = -111,
    TYPH_ERR_MISSING_ARG   = -112,
    TYPH_ERR_INVALID_OP    = -113,
    TYPH_ERR_UNKNOWN_MODE  = -114
};



/** Specify a send or receive operation. */
enum TYPH_Sendrecv : int {
    TYPH_SENDRECV_SEND = 1001,
    TYPH_SENDRECV_RECV = 1002
};



/**
 * Each send/recv op for each ghost layer in each quant requires a key, which
 * defines the memory displacements necessary for creating an MPI derived type.
 */
struct Key
{
    static constexpr int SERIALISATION_MARKER = 0x0;

    int layer = -1;         //!< ghost layer that key is concerned with
    int proc = -1;          //!< proc to send to or receive from

    int *list = nullptr;    //!< list of elements or nodes (indices within layer)
    int list_len = 0;       //!< size of element/node list

    //! Optional block lengths if variable amounts of data to be communicated
    //! per cell/node
    int *block_lens = nullptr;

    Key *parent = nullptr;  //!< key this key is derived from (for auxiliary
                            //!  data)

    Key *next   = nullptr;  //!< next key in key set linked list
    Key *prev   = nullptr;  //!< previous key in key set linked list
};



/**
 * The collection of keys necessary for communicating a quant form a key set. If
 * different quants have the same structure in memory, they may share a key set.
 */
struct Key_Set
{
    static constexpr int SERIALISATION_MARKER = 0x1;

    //std::string name;           // key set name

    TYPH_Centring centring;     // Type of variable centring (NODE/CELL)
    TYPH_Auxiliary aux;         // Type of auxiliary data

    int stride = 0;             // # parts object is divided into (for 2D etc. arrays)
    int layer_min = -1;         // smallest ghost layer index
    int layer_max = -1;         // largest ghost layer index
    int partition_id = -1;      // ID of partition to use with keys

    int num_send = 0;           // # send keys in set
    int num_recv = 0;           // # recv keys in set

    Key *send_keys = nullptr;   // linked list of send keys
    Key *recv_keys = nullptr;   // linked list of recv keys
};



/** TODO */
struct Quant_Info
{
    static constexpr int SERIALISATION_MARKER = 0x4;

    int quant_id = -1;      // Quant ID
    int recv_quant_id = -1;
    int key_set_id = -1;    // Associated key set ID

    int layer_min = -1;     // Default ghost layer range
    int layer_max = -1;     // Default ghost layer range

    int quant_size = 0;
    int nrepeat = 0;
    int stride = 0;

    MPI_Datatype *old_mpi_type = nullptr;
    MPI_Datatype *new_mpi_type = nullptr;

    Key_Set    *key_set = nullptr;   // Key set for this quant
    Quant_Info *next    = nullptr;
};



/** TODO */
struct Schedule_Part
{
    static constexpr int SERIALISATION_MARKER = 0x5;

    MPI_Datatype new_mpi_type = MPI_DATATYPE_NULL;
    int quant_size = 0;
    int nrepeat = 0;
    int stride = 0;
    Key *key = nullptr;             // offset array and length
    MPI_Datatype *old_mpi_type = nullptr;

    MPI_Aint *address = nullptr;   // Start address of item
};



/** TODO */
struct Schedule
{
    static constexpr int SERIALISATION_MARKER = 0x3;

    static int max_parts;           //!< Initialised to 0 in schedule.cpp
    static int max_procs;           //!< Initialised to 0 in schedule.cpp

    int num_send = 0;               //!< \# procs to send to
    int *send_procs;                //!< Procs to send to
    MPI_Datatype *mpi_send_types;   //!< MPI types of send data
    MPI_Request *send_requests;     //!< MPI send requests
    int *send_nparts;               //!< \# parts to send to each proc
    int *send_start;                //!< Start index in parts array for send to i

    int num_recv = 0;               //!< \# procs to recv from
    int *recv_procs;                //!< Procs to recv from
    MPI_Datatype *mpi_recv_types;   //!< MPI types of recv data
    MPI_Request *recv_requests;     //!< MPI recv requests
    int *recv_nparts;               //!< \# parts to recv from each proc
    int *recv_start;                //!< Start index in parts array for recv from i

    Schedule_Part *parts;           //!< Array of all parts to/from all procs for
                                    //!  this proc
};



/** TODO */
struct Phase
{
    /** Mark a phase in serialised output. */
    static constexpr int SERIALISATION_MARKER = 0x2;

    std::string name;                   //!< Phase name

    bool is_pure = true;                //!< Is this a pure or auxiliary phase

    int *quant_id_list = nullptr;       //!< List of quant IDs
    int num_quants = 0;                 //!< \# quants in the phase

    int key_set_id = -1;                //!< Key set ID

    int num_layers = 0;                 //!< \# ghost layers
    int layer_min  = 0;                 //!< Minimum ghost layer (inclusive)
    int layer_max  = 0;                 //!< Maximum ghost layer (inclusive)

    Quant_Info *quant_info = nullptr;   //!< Linked list of quant metadata
    Schedule   *schedule   = nullptr;   //!< Comms schedule

    bool is_built     = false;          //!< Schedule built?
    bool is_committed = false;          //!< Phase MPI committed?
};



/** TODO */
struct Quant
{
    /** Mark a quant in serialised output. */
    static constexpr int SERIALISATION_MARKER = 0x6;

    std::string name;

    int quant_data_id = -1;
    int num_layers = 0;
    TYPH_Centring centring;
    TYPH_Datatype datatype;
    TYPH_Auxiliary aux;         // Handle for type of auxiliary data
    bool is_pure = true;
    bool is_aux = false;

    int rank = -1;              // Rank of variable
    int mesh_dim = -1;          // Mesh-based dimension
    int stride = 0;             // Multi-dimensional array stride
    int *dims = nullptr;        // Non-mesh-based dimensions if rank > 1

    // Array index bounds. These are just set to [0,dim) at the moment.
    int *lower_bound = nullptr;
    int *upper_bound = nullptr;

    MPI_Aint quant_address = 0; // Start address of quant
    MPI_Datatype mpi_datatype = MPI_DATATYPE_NULL;
};

/** @} */



// -----------------------------------------------------------------------------
// Internal constants
// -----------------------------------------------------------------------------
/**
 * \addtogroup internal
 *
 * @{
 */

// Error classes
int constexpr ERR_USER = 101;   //!< User error
int constexpr ERR_MEM  = 102;   //!< Memory error
int constexpr ERR_MPI  = 103;   //!< MPI error
int constexpr ERR_INT  = 104;   //!< Internal error
int constexpr ERR_APP  = 105;   //!< For when application calls abort

MPI_Aint constexpr TYPH_NULL_ADDR = 0;  //!< Null quantity address

/** @} */

} // namespace _TYPH_Internal



#endif // TYPHON_TYPES_H
