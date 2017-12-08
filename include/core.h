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
#ifndef TYPHON_CORE_H
#define TYPHON_CORE_H

#include <functional>
#include <vector>



namespace _TYPH_Internal {

/**
 * \addtogroup internal
 *
 * @{
 */

/** Function pointer to execute on shutdown.*/
using Kill_Func = std::function<int()>;



/** Stores basic information about the underlying MPI configuration. */
struct MPI_Runtime
{
    MPI_Comm comm   = MPI_COMM_NULL;    //!< MPI communicator
    MPI_Info info   = MPI_INFO_NULL;    //!< MPI info

    int size        = 0;                //!< Communicator size
    int rank        = MPI_PROC_NULL;    //!< Communicator rank
    int min_rank    = MPI_PROC_NULL;    //!< Smallest rank in communicator
    int max_rank    = MPI_PROC_NULL;    //!< Largest rank in communicator
    int master_rank = MPI_PROC_NULL;    //!< Master rank in communicator
    int error       = MPI_SUCCESS;      //!< MPI error

    bool is_master   = false;           //!< Is this rank master?
    bool initialised = false;           //!< Has MPI been initialised?
    bool finalised   = false;           //!< Has MPI been finalised?
};



/** Singleton keeping track of Typhon runtime information. */
class Core
{
public:
    /** Initialise the global singleton instance. */
    static void Init_Singleton();

    /** Get a pointer to the global singleton instances. */
    static Core *Get_Singleton() { return singleton; }

    /** Destruct the global singleton instance. */
    static void Kill_Singleton();

    /** Return a pointer to the underlying MPI_Runtime. */
    MPI_Runtime *Get_MPI_Runtime() { return mpi_runtime; }

    /** Is Typhon initialised? */
    bool Is_Initialised() const
        { return mpi_runtime != nullptr && mpi_runtime->initialised; }

    /** Initialise Typhon. */
    int Initialise(MPI_Comm *comm);

    /** Shutdown Typhon. */
    int Kill(bool finalise);

    /** Add a hook to be executed on shutdown. */
    int Kill_Insert(Kill_Func kill_func);

    /** Synchronise processors. */
    int Comms_Barrier();

    /** Abort execution. */
    int Comms_Abort(int abort_code);

private:
    static Core *singleton; //!< Global singleton instance

    /** Constructor. Construction only allowed through Get_Singleton(). */
    explicit Core();

    /** Destructor. Destruction only allowed through Kill_Singleton(). */
    ~Core();

    /** Unimplemented copy-constructor---disallow copying. */
    Core(Core const &rhs);

    /** Unimplemented operator=---disallow copying. */
    Core &operator=(Core const &rhs);

    MPI_Runtime *mpi_runtime;           //!< MPI runtime information
    std::vector<Kill_Func> kill_store;  //!< Functions to call on shutdown

    int Comms_Init(MPI_Comm *comm);
    int Comms_Kill(bool finalise);
    int Init_MPI_Runtime(MPI_Comm *comm);
};

#define TYPH_CORE _TYPH_Internal::Core::Get_Singleton()

/** @} */

} // namespace _TYPH_Internal

// -----------------------------------------------------------------------------
// Public Typhon API - core functions
// -----------------------------------------------------------------------------
/**
 * \addtogroup typhon
 *
 * @{
 */

/** Initialise Typhon. */
int
TYPH_Init(
        MPI_Comm *comm = nullptr);

/** Shutdown Typhon. */
int
TYPH_Kill(
        bool finalise = true);

/** Abort Typhon in the case of an error. */
int
TYPH_Abort(
        int abort_code);

/** Get the current processor count. */
int
TYPH_Get_Size(
        int *size);

/** Get the rank of the calling processor. */
int
TYPH_Get_Rank(
        int *rank);

/** Return whether or not the calling processor is the master rank. */
bool
TYPH_Is_Master();

/** Synchronise all processors. */
int
TYPH_Barrier();

/** Return timestamp. */
double
TYPH_Get_Time();

/** Set the argument to the current MPI communicator. */
int
TYPH_Set_Comm(
        MPI_Comm *comm);

/** Set the argument to the MPI self communicator. */
int
TYPH_Set_Comm_Self(
        MPI_Comm *comm);

/** @} */



#endif // TYPHON_CORE_H
