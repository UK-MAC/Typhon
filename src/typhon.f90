module typhon

use, intrinsic :: iso_fortran_env, ONLY: real64, int32
use, intrinsic :: iso_c_binding

use mpi

implicit none

!
! Kinds for Fortran convenience wrappers
!

integer, parameter :: f64 = real64
integer, parameter :: i32 = int32
integer, parameter :: z32 = int32


!
! Typhon public API constants
!

integer(i32), parameter :: typh_mesh_dim = -88

integer(i32), parameter :: typh_success = 0
integer(i32), parameter :: typh_fail    = -1

integer(i32), parameter :: typh_pure = 1

enum, bind(c)
    enumerator :: typh_datatype_null = 0, &
                  typh_datatype_real = 1, &
                  typh_datatype_integer = 2, &
                  typh_datatype_logical = 3
end enum

enum, bind(c)
    enumerator :: typh_auxiliary_none = -9099
end enum

enum, bind(c)
    enumerator :: typh_ghosts_zero  = 0, &
                  typh_ghosts_one   = 1, &
                  typh_ghosts_two   = 2, &
                  typh_ghosts_three = 3
end enum

enum, bind(c)
    enumerator :: typh_shape_quad = 4
end enum

enum, bind(c)
    enumerator :: typh_centring_node = 2001, &
                  typh_centring_cell = 2002
end enum

enum, bind(c)
    enumerator :: typh_keytype_cell        = 1, &
                  typh_keytype_node        = 2, &
                  typh_keytype_cell_corner = 3
end enum

enum, bind(c)
    enumerator :: typh_op_sum  = 1001, &
                  typh_op_prod = 1002, &
                  typh_op_max  = 1003, &
                  typh_op_min  = 1004, &
                  typh_op_or   = 1011, &
                  typh_op_xor  = 1012, &
                  typh_op_and  = 1013
end enum


!
! Bind to Typhon public API functions
!

interface
    function typh_free_c(ptr) bind(c, name="TYPH_Free")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        integer(c_int) :: typh_free_c
        type(c_ptr), value :: ptr
    end function typh_free_c

    function typh_init_c(comm) bind(c, name="TYPH_Init")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        integer(c_int) :: typh_init_c
        type(c_ptr), value :: comm
    end function typh_init_c

    function typh_init_f_c(comm) bind(c, name="TYPH_Init_F")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_init_f_c
        integer :: comm
    end function typh_init_f_c

    function typh_kill_c(finalise) bind(c, name="TYPH_Kill")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_kill_c
        integer(c_int), value :: finalise
    end function typh_kill_c

    function typh_abort_c(abort_code) bind(c, name="TYPH_Abort")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_abort_c
        integer(c_int), value :: abort_code
    end function typh_abort_c

    function typh_get_size_c(size) bind(c, name="TYPH_Get_Size")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_get_size_c
        integer(c_int) :: size
    end function typh_get_size_c

    function typh_get_rank_c(rank) bind(c, name="TYPH_Get_Rank")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_get_rank_c
        integer(c_int) :: rank
    end function typh_get_rank_c

    function typh_is_master_c() bind(c, name="TYPH_Is_Master")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_is_master_c
    end function typh_is_master_c

    function typh_barrier_c() bind(c, name="TYPH_Barrier")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_barrier_c
    end function typh_barrier_c

    function typh_get_time_c() bind(c, name="TYPH_Get_Time")
        use, intrinsic :: iso_c_binding, only: c_double
        real(c_double) :: typh_get_time_c
    end function typh_get_time_c

    function typh_set_comm_c(comm) bind(c, name="TYPH_Set_Comm_F")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_set_comm_c
        integer :: comm
    end function typh_set_comm_c

    function typh_set_comm_self_c(comm) bind(c, name="TYPH_Set_Comm_Self_F")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_set_comm_self_c
        integer :: comm
    end function typh_set_comm_self_c

    function typh_start_register_c() bind(c, name="TYPH_Start_Register")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_start_register_c
    end function typh_start_register_c

    function typh_finish_register_c() bind(c, name="TYPH_Finish_Register")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_finish_register_c
    end function typh_finish_register_c

    function typh_is_registering_c(is_registering) &
            bind(c, name="TYPH_Is_Registering")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_is_registering_c
        integer(c_int) :: is_registering
    end function typh_is_registering_c

    function typh_add_phase_c(phase_id, name, nghosts, is_pure, key_set_id, &
            ghosts_min) &
            bind(c, name="TYPH_Add_Phase")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_add_phase_c

        integer(c_int) :: phase_id
        type(c_ptr), value :: name
        integer(kind(typh_ghosts_zero)), value :: nghosts
        integer(c_int), value :: is_pure
        integer(c_int), value :: key_set_id
        integer(c_int), value :: ghosts_min

    end function typh_add_phase_c

    function typh_add_quant_c(quant_id, name, num_ghosts, datatype, &
            centring, is_pure, aux, dims, num_dims) &
            bind(c, name="TYPH_Add_Quant")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_add_quant_c

        integer(c_int) :: quant_id
        type(c_ptr), value :: name
        integer(kind(typh_ghosts_zero)), value :: num_ghosts
        integer(kind(typh_datatype_null)), value :: datatype
        integer(kind(typh_centring_cell)), value :: centring
        integer(c_int), value :: is_pure
        integer(kind(typh_auxiliary_none)), value :: aux
        type(c_ptr), value :: dims
        integer(c_int), value :: num_dims

    end function typh_add_quant_c

    function typh_distribute_mesh_c(nel_local, nel_global, nnd_global, &
            nd_per_el, num_layers, conn_data, partition, layer_nel, &
            layer_nnd, total_nel, total_nnd, el_loc_glob, nd_loc_glob, &
            el_region, el_material, el_nd, el_owner, nd_owner) &
            bind(c, name="TYPH_Distribute_Mesh")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        integer(c_int) :: typh_distribute_mesh_c
        integer(c_int), value :: nel_local
        integer(c_int), value :: nel_global
        integer(c_int), value :: nnd_global
        integer(c_int), value :: nd_per_el
        integer(c_int), value :: num_layers
        type(c_ptr), value :: conn_data
        type(c_ptr), value :: partition
        type(c_ptr), value :: layer_nel
        type(c_ptr), value :: layer_nnd
        integer(c_int) :: total_nel
        integer(c_int) :: total_nnd
        type(c_ptr) :: el_loc_glob
        type(c_ptr) :: nd_loc_glob
        type(c_ptr) :: el_region
        type(c_ptr) :: el_material
        type(c_ptr) :: el_nd
        type(c_ptr) :: el_owner
        type(c_ptr) :: nd_owner
    end function typh_distribute_mesh_c

    function typh_set_partition_info_c(partition_id, el_shape, num_layers, &
            num_el_total, num_nd_total, el_owner, nd_owner, el_loc_glob, &
            nd_loc_glob, connectivity, name) &
            bind(c, name="TYPH_Set_Partition_Info")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        integer(c_int) :: typh_set_partition_info_c
        integer(c_int) :: partition_id
        integer(kind(typh_shape_quad)), value :: el_shape
        integer(c_int), value :: num_layers
        type(c_ptr), value :: num_el_total
        type(c_ptr), value :: num_nd_total
        type(c_ptr), value :: el_owner
        type(c_ptr), value :: nd_owner
        type(c_ptr), value :: el_loc_glob
        type(c_ptr), value :: nd_loc_glob
        type(c_ptr), value :: connectivity
        type(c_ptr), value :: name
    end function typh_set_partition_info_c

    function typh_create_key_set_c(key_type, layer_min, layer_max, &
            partition_id, key_set_id) &
            bind(c, name="TYPH_Create_Key_Set")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr
        integer(c_int) :: typh_create_key_set_c
        integer(kind(typh_keytype_cell)), value :: key_type
        integer(c_int), value :: layer_min
        integer(c_int), value :: layer_max
        integer(c_int), value :: partition_id
        integer(c_int) :: key_set_id
    end function typh_create_key_set_c

    function typh_add_quant_to_phase_c(phase_id, quant_id, recv_quant_id, &
            key_set_id, ghosts_min, ghosts_max) &
            bind(c, name="TYPH_Add_Quant_To_Phase")
        use, intrinsic :: iso_c_binding, only: c_int

        integer(c_int) :: typh_add_quant_to_phase_c

        integer(c_int), value :: phase_id
        integer(c_int), value :: quant_id
        integer(c_int), value :: recv_quant_id
        integer(c_int), value :: key_set_id
        integer(c_int), value :: ghosts_min
        integer(c_int), value :: ghosts_max

    end function typh_add_quant_to_phase_c

    function typh_set_quant_address_c(quant_id, data, dims, rank) &
            bind(c, name="TYPH_Set_Quant_Address")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_set_quant_address_c

        integer(c_int), value :: quant_id
        type(c_ptr), value :: data
        type(c_ptr), value :: dims
        integer(c_int), value :: rank

    end function typh_set_quant_address_c

    function typh_exchange_c(phase_id) &
            bind(c, name="TYPH_Exchange")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_exchange_c
        integer(c_int), value :: phase_id
    end function typh_exchange_c

    function typh_start_exchange_c(phase_id) &
            bind(c, name="TYPH_Start_Exchange")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_start_exchange_c
        integer(c_int), value :: phase_id
    end function typh_start_exchange_c

    function typh_finish_exchange_c(phase_id) &
            bind(c, name="TYPH_Finish_Exchange")
        use, intrinsic :: iso_c_binding, only: c_int
        integer(c_int) :: typh_finish_exchange_c
        integer(c_int), value :: phase_id
    end function typh_finish_exchange_c

    function typh_gather_i_c(in, dims, rank, out) &
            bind(c, name="TYPH_Gather_i")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_gather_i_c

        type(c_ptr), value :: in
        type(c_ptr), value :: dims
        integer(c_int), value :: rank
        type(c_ptr), value :: out

    end function typh_gather_i_c

    function typh_gather_d_c(in, dims, rank, out) &
            bind(c, name="TYPH_Gather_d")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_gather_d_c

        type(c_ptr), value :: in
        type(c_ptr), value :: dims
        integer(c_int), value :: rank
        type(c_ptr), value :: out

    end function typh_gather_d_c

    function typh_gather_z_c(in, dims, rank, out) &
            bind(c, name="TYPH_Gather_z")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_gather_z_c

        type(c_ptr), value :: in
        type(c_ptr), value :: dims
        integer(c_int), value :: rank
        type(c_ptr), value :: out

    end function typh_gather_z_c

    function typh_gather_c(datatype, in, dims, rank, out) &
            bind(c, name="TYPH_Gather")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_gather_c

        integer(kind(typh_datatype_null)), value :: datatype
        type(c_ptr), value :: in
        type(c_ptr), value :: dims
        integer(c_int), value :: rank
        type(c_ptr), value :: out

    end function typh_gather_c

    function typh_reduce_i_c(in, dims, rank, out, op) &
            bind(c, name="TYPH_Reduce_i")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_reduce_i_c

        type(c_ptr), value :: in
        type(c_ptr), value :: dims
        integer(c_int), value :: rank
        type(c_ptr), value :: out
        integer(kind(typh_op_sum)), value :: op

    end function typh_reduce_i_c

    function typh_reduce_d_c(in, dims, rank, out, op) &
            bind(c, name="TYPH_Reduce_d")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_reduce_d_c

        type(c_ptr), value :: in
        type(c_ptr), value :: dims
        integer(c_int), value :: rank
        type(c_ptr), value :: out
        integer(kind(typh_op_sum)), value :: op

    end function typh_reduce_d_c

    function typh_reduce_z_c(in, dims, rank, out, op) &
            bind(c, name="TYPH_Reduce_z")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_reduce_z_c

        type(c_ptr), value :: in
        type(c_ptr), value :: dims
        integer(c_int), value :: rank
        type(c_ptr), value :: out
        integer(kind(typh_op_sum)), value :: op

    end function typh_reduce_z_c

    function typh_reduce_c(datatype, in, dims, rank, out, op) &
            bind(c, name="TYPH_Reduce")
        use, intrinsic :: iso_c_binding, only: c_int, c_ptr

        integer(c_int) :: typh_reduce_c

        integer(kind(typh_datatype_null)), value :: datatype
        type(c_ptr), value :: in
        type(c_ptr), value :: dims
        integer(c_int), value :: rank
        type(c_ptr), value :: out
        integer(kind(typh_op_sum)), value :: op

    end function typh_reduce_c
end interface

interface typh_init
    module procedure typh_init_bare, typh_init_comm
end interface typh_init

private

public f64, i32, z32, typh_mesh_dim, typh_success, typh_fail, typh_pure, &
    typh_datatype_null, typh_datatype_real, typh_datatype_integer, &
    typh_datatype_logical, &
    typh_auxiliary_none, &
    typh_ghosts_zero, typh_ghosts_one, typh_ghosts_two, typh_ghosts_three, &
    typh_shape_quad, &
    typh_centring_node, typh_centring_cell, &
    typh_keytype_cell, typh_keytype_node, typh_keytype_cell_corner, &
    typh_free, typh_init, typh_kill, typh_abort, typh_get_size, typh_get_rank, &
    typh_is_master, typh_barrier, typh_get_time, &
    typh_set_comm, typh_set_comm_self, &
    typh_start_register, typh_finish_register, typh_is_registering, &
    typh_add_phase, typh_add_quant, &
    typh_distribute_mesh, typh_set_partition_info, &
    typh_create_key_set, typh_add_quant_to_phase, typh_set_quant_address, &
    typh_exchange, typh_start_exchange, typh_finish_exchange, &
    typh_gather_i, typh_gather_d, typh_gather_z, typh_gather, &
    typh_reduce_i, typh_reduce_d, typh_reduce_z, typh_reduce


contains

!
! Provide Fortran convenience wrappers for the above bindings
!

subroutine typh_free(ptr, typh_err)

    type(c_ptr), intent(in) :: ptr
    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_free_c(ptr), kind=i32)

end subroutine typh_free


subroutine typh_init_bare(typh_err)

    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_init_c(c_null_ptr), kind=i32)

end subroutine typh_init_bare


subroutine typh_init_comm(comm, typh_err)

    integer, intent(inout) :: comm
    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_init_f_c(comm), kind=i32)

end subroutine typh_init_comm


subroutine typh_kill(finalise, typh_err)

    logical(z32), intent(in) :: finalise
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: finalise_c
    finalise_c = merge(1_c_int, 0_c_int, finalise)

    typh_err = int(typh_kill_c(finalise_c), kind=i32)

end subroutine typh_kill


subroutine typh_abort(abort_code, typh_err)

    integer(i32), intent(in) :: abort_code
    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_abort_c(int(abort_code, kind=c_int)), kind=i32)

end subroutine typh_abort


subroutine typh_get_size(size, typh_err)

    integer(i32), intent(out) :: size
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: size_c
    typh_err = int(typh_get_size_c(size_c), kind=i32)
    size = int(size_c, kind=i32)

end subroutine typh_get_size


subroutine typh_get_rank(rank, typh_err)

    integer(i32), intent(out) :: rank
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    typh_err = int(typh_get_rank_c(rank_c), kind=i32)
    rank = int(rank_c, kind=i32)

end subroutine typh_get_rank


subroutine typh_is_master(is_master, typh_err)

    logical(z32), intent(out) :: is_master
    integer(i32), intent(out) :: typh_err

    is_master = (typh_is_master_c() == 1_c_int)
    typh_err = typh_success

end subroutine typh_is_master


subroutine typh_barrier(typh_err)

    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_barrier_c(), kind=i32)

end subroutine typh_barrier


subroutine typh_get_time(time, typh_err)

    real(f64), intent(out) :: time
    integer(i32), intent(out) :: typh_err

    time = real(typh_get_time_c(), kind=f64)
    typh_err = typh_success

end subroutine typh_get_time


subroutine typh_set_comm(comm, typh_err)

    integer, intent(inout) :: comm
    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_set_comm_c(comm), kind=i32)

end subroutine typh_set_comm


subroutine typh_set_comm_self(comm, typh_err)

    integer, intent(inout) :: comm
    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_set_comm_self_c(comm), kind=i32)

end subroutine typh_set_comm_self


subroutine typh_start_register(typh_err)

    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_start_register_c(), kind=i32)

end subroutine typh_start_register


subroutine typh_finish_register(typh_err)

    integer(i32), intent(out) :: typh_err

    typh_err = int(typh_finish_register_c(), kind=i32)

end subroutine typh_finish_register


subroutine typh_is_registering(is_registering, typh_err)

    integer(i32), intent(out) :: is_registering
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: is_registering_c
    typh_err = int(typh_is_registering_c(is_registering_c), kind=i32)
    is_registering = int(is_registering_c, kind=i32)

end subroutine typh_is_registering


subroutine typh_add_phase(phase_id, name, nghosts, is_pure, key_set_id, &
        ghosts_min, typh_err)

    integer(i32), intent(out) :: phase_id
    character(len=*), intent(in) :: name
    integer(kind(typh_ghosts_zero)), intent(in) :: nghosts
    logical(z32), intent(in) :: is_pure
    integer(i32), intent(in) :: key_set_id
    integer(i32), intent(in) :: ghosts_min
    integer(i32), intent(out) :: typh_err

    character(len=len(name)+1), target :: tmp_name
    type(c_ptr) :: name_c
    integer(c_int) :: phase_id_c, is_pure_c, key_set_id_c, &
        ghosts_min_c

    tmp_name = name//c_null_char

    name_c = c_loc(tmp_name)
    is_pure_c = merge(1_c_int, 0_c_int, is_pure)
    key_set_id_c = int(key_set_id, kind=c_int)
    ghosts_min_c = int(ghosts_min, kind=c_int)

    typh_err = int(typh_add_phase_c(phase_id_c, name_c, nghosts, is_pure_c, &
        key_set_id_c, ghosts_min_c), kind=c_int)
    phase_id = int(phase_id_c, kind=i32)

end subroutine typh_add_phase


subroutine typh_add_quant(quant_id, name, num_ghosts, datatype, &
        centring, is_pure, aux, dims, num_dims, typh_err)

    integer(i32), intent(out) :: quant_id
    character(len=*), intent(in) :: name
    integer(kind(typh_ghosts_zero)), intent(in) :: num_ghosts
    integer(kind(typh_datatype_null)), intent(in) :: datatype
    integer(kind(typh_centring_cell)), intent(in) :: centring
    logical(z32), intent(in) :: is_pure
    integer(kind(typh_auxiliary_none)), intent(in) :: aux
    integer(i32), intent(in) :: num_dims
    integer(i32), dimension(num_dims), target, intent(in) :: dims
    integer(i32), intent(out) :: typh_err

    character(len=len(name)+1), target :: tmp_name
    integer(c_int) :: quant_id_c, is_pure_c, num_dims_c
    type(c_ptr) :: name_c, dims_c

    tmp_name = name//c_null_char

    is_pure_c = merge(1_c_int, 0_c_int, is_pure)
    num_dims_c = int(num_dims, kind=c_int)
    name_c = c_loc(tmp_name)
    dims_c = c_loc(dims)

    typh_err = int(typh_add_quant_c(quant_id_c, name_c, num_ghosts, datatype, &
        centring, is_pure_c, aux, dims_c, num_dims_c), kind=i32)

    quant_id = int(quant_id_c, kind=i32)

end subroutine typh_add_quant


subroutine typh_distribute_mesh(nel_local, nnd_local, nel_global, nnd_global, &
        nd_per_el, num_layers, conn_data, partition, layer_nel, layer_nnd, &
        total_nel, total_nnd, el_loc_glob, nd_loc_glob, el_region, &
        el_material, el_nd, el_owner, nd_owner, typh_err)

    integer(i32), intent(in) :: nel_local
    integer(i32), intent(in) :: nnd_local
    integer(i32), intent(in) :: nel_global
    integer(i32), intent(in) :: nnd_global
    integer(i32), intent(in) :: nd_per_el
    integer(i32), intent(in) :: num_layers
    integer(i32), dimension(3+nd_per_el,nel_local), target, intent(in) :: &
        conn_data
    integer(i32), dimension(nel_local), target, intent(in) :: partition
    integer(i32), dimension(num_layers+1), target, intent(out) :: layer_nel
    integer(i32), dimension(num_layers+1), target, intent(out) :: layer_nnd
    integer(i32), intent(out) :: total_nel
    integer(i32), intent(out) :: total_nnd
    integer(i32), pointer, intent(out) :: el_loc_glob(:)
    integer(i32), pointer, intent(out) :: nd_loc_glob(:)
    integer(i32), pointer, intent(out) :: el_region(:)
    integer(i32), pointer, intent(out) :: el_material(:)
    integer(i32), pointer, intent(out) :: el_nd(:,:)
    integer(i32), pointer, intent(out) :: el_owner(:,:)
    integer(i32), pointer, intent(out) :: nd_owner(:,:)
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: nel_local_c, nel_global_c, nnd_global_c, nd_per_el_c, &
        num_layers_c, total_nel_c, total_nnd_c
    type(c_ptr) :: conn_data_c, partition_c, layer_nel_c, layer_nnd_c, &
        el_loc_glob_c, nd_loc_glob_c, el_region_c, el_material_c, &
        el_nd_c, el_owner_c, nd_owner_c

    nel_local_c = int(nel_local, kind=c_int)
    nel_global_c = int(nel_global, kind=c_int)
    nnd_global_c = int(nnd_global, kind=c_int)
    nd_per_el_c = int(nd_per_el, kind=c_int)
    num_layers_c = int(num_layers, kind=c_int)
    conn_data_c = c_loc(conn_data)
    partition_c = c_loc(partition)
    layer_nel_c = c_loc(layer_nel)
    layer_nnd_c = c_loc(layer_nnd)

    typh_err = int(typh_distribute_mesh_c(nel_local_c, nel_global_c, &
        nnd_global_c, nd_per_el_c, num_layers_c, conn_data_c, &
        partition_c, layer_nel_c, layer_nnd_c, total_nel_c, total_nnd_c, &
        el_loc_glob_c, nd_loc_glob_c, el_region_c, el_material_c, el_nd_c, &
        el_owner_c, nd_owner_c), kind=i32)

    total_nel = int(total_nel_c, kind=i32)
    total_nnd = int(total_nnd_c, kind=i32)

    call c_f_pointer(el_loc_glob_c, el_loc_glob, [nel_local])
    call c_f_pointer(nd_loc_glob_c, nd_loc_glob, [nnd_local])
    call c_f_pointer(el_region_c, el_region, [nel_local])
    call c_f_pointer(el_material_c, el_material, [nel_local])
    call c_f_pointer(el_nd_c, el_nd, [4, nel_local])
    call c_f_pointer(el_owner_c, el_owner, [2, nel_local])
    call c_f_pointer(nd_owner_c, nd_owner, [2, nnd_local])

end subroutine typh_distribute_mesh


subroutine typh_set_partition_info(partition_id, nel_local, nnd_local, &
        el_shape, num_layers, num_el_total, num_nd_total, el_owner, &
        nd_owner, el_loc_glob, nd_loc_glob, el_nd, name, typh_err)

    integer(i32), intent(out) :: partition_id
    integer(i32), intent(in) :: nel_local
    integer(i32), intent(in) :: nnd_local
    integer(kind(typh_ghosts_zero)), intent(in) :: el_shape
    integer(i32), intent(in) :: num_layers
    integer(i32), dimension(num_layers+1), target, intent(in) :: num_el_total
    integer(i32), dimension(num_layers+1), target, intent(in) :: num_nd_total
    integer(i32), dimension(2,nel_local), target, intent(in) :: el_owner
    integer(i32), dimension(2,nnd_local), target, intent(in) :: nd_owner
    integer(i32), dimension(nel_local), target, intent(in) :: el_loc_glob
    integer(i32), dimension(nnd_local), target, intent(in) :: nd_loc_glob
    integer(i32), dimension(4,nel_local), target, intent(in) :: el_nd
    character(len=*), intent(in) :: name
    integer(i32), intent(out) :: typh_err

    character(len=len(name)+1), target :: tmp_name
    integer(c_int) :: partition_id_c, el_shape_c, num_layers_c
    type(c_ptr) :: num_el_total_c, num_nd_total_c, el_owner_c, nd_owner_c, &
        el_loc_glob_c, nd_loc_glob_c, el_nd_c, name_c

    tmp_name = name//c_null_char

    el_shape_c = int(el_shape, kind=c_int)
    num_layers_c = int(num_layers, kind=c_int)

    num_el_total_c = c_loc(num_el_total)
    num_nd_total_c = c_loc(num_nd_total)
    el_owner_c = c_loc(el_owner)
    nd_owner_c = c_loc(nd_owner)
    el_loc_glob_c = c_loc(el_loc_glob)
    nd_loc_glob_c = c_loc(nd_loc_glob)
    el_nd_c = c_loc(el_nd)
    name_c = c_loc(tmp_name)

    typh_err = int(typh_set_partition_info_c(partition_id_c, el_shape_c, &
        num_layers_c, num_el_total_c, num_nd_total_c, el_owner_c, nd_owner_c, &
        el_loc_glob_c, nd_loc_glob_c, el_nd_c, name_c), kind=i32)

    partition_id = int(partition_id_c, kind=i32)

end subroutine typh_set_partition_info


subroutine typh_create_key_set(key_type, layer_min, layer_max, partition_id, &
        key_set_id, typh_err)

    integer(kind(typh_keytype_node)), intent(in) :: key_type
    integer(i32), intent(in) :: layer_min
    integer(i32), intent(in) :: layer_max
    integer(i32), intent(in) :: partition_id
    integer(i32), intent(out) :: key_set_id
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: layer_min_c, layer_max_c, partition_id_c, key_set_id_c

    layer_min_c = int(layer_min, kind=c_int)
    layer_max_c = int(layer_max, kind=c_int)
    partition_id_c = int(partition_id, kind=c_int)

    typh_err = int(typh_create_key_set_c(key_type, layer_min_c, layer_max_c, &
        partition_id_c, key_set_id_c), kind=i32)

    key_set_id = int(key_set_id_c, kind=i32)

end subroutine typh_create_key_set


subroutine typh_add_quant_to_phase(phase_id, quant_id, recv_quant_id, &
        key_set_id, ghosts_min, ghosts_max, typh_err)

    integer(i32), intent(in) :: phase_id
    integer(i32), intent(in) :: quant_id
    integer(i32), intent(in) :: recv_quant_id
    integer(i32), intent(in) :: key_set_id
    integer(i32), intent(in) :: ghosts_min
    integer(i32), intent(in) :: ghosts_max
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: phase_id_c, quant_id_c, recv_quant_id_c, key_set_id_c, &
        ghosts_min_c, ghosts_max_c

    phase_id_c = int(phase_id, kind=c_int)
    quant_id_c = int(quant_id, kind=c_int)
    recv_quant_id_c = int(recv_quant_id, kind=c_int)
    key_set_id_c = int(key_set_id, kind=c_int)
    ghosts_min_c = int(ghosts_min, kind=c_int)
    ghosts_max_c = int(ghosts_max, kind=c_int)

    typh_err = int(typh_add_quant_to_phase_c(phase_id_c, quant_id_c, &
        recv_quant_id_c, key_set_id_c, ghosts_min_c, ghosts_max_c), kind=i32)

end subroutine typh_add_quant_to_phase


subroutine typh_set_quant_address(quant_id, data, dims, rank, typh_err)

    integer(i32), intent(in) :: quant_id
    type(c_ptr), intent(in) :: data
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: quant_id_c, rank_c
    type(c_ptr) :: dims_c

    quant_id_c = int(quant_id, kind=c_int)
    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_set_quant_address_c(quant_id_c, data, dims_c, &
        rank_c), kind=i32)

end subroutine typh_set_quant_address


subroutine typh_exchange(phase_id, typh_err)

    integer(i32), intent(in) :: phase_id
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: phase_id_c

    phase_id_c = int(phase_id, kind=c_int)

    typh_err = int(typh_exchange_c(phase_id_c), kind=i32)

end subroutine typh_exchange


subroutine typh_start_exchange(phase_id, typh_err)

    integer(i32), intent(in) :: phase_id
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: phase_id_c

    phase_id_c = int(phase_id, kind=c_int)

    typh_err = int(typh_start_exchange_c(phase_id_c), kind=i32)

end subroutine typh_start_exchange


subroutine typh_finish_exchange(phase_id, typh_err)

    integer(i32), intent(in) :: phase_id
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: phase_id_c

    phase_id_c = int(phase_id, kind=c_int)

    typh_err = int(typh_finish_exchange_c(phase_id_c), kind=i32)

end subroutine typh_finish_exchange


subroutine typh_gather_i(in, dims, rank, out, typh_err)

    type(c_ptr), intent(in) :: in
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    type(c_ptr), intent(out) :: out
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    type(c_ptr) :: dims_c

    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_gather_i_c(in, dims_c, rank_c, out), kind=i32)

end subroutine typh_gather_i


subroutine typh_gather_d(in, dims, rank, out, typh_err)

    type(c_ptr), intent(in) :: in
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    type(c_ptr), intent(out) :: out
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    type(c_ptr) :: dims_c

    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_gather_d_c(in, dims_c, rank_c, out), kind=i32)

end subroutine typh_gather_d


subroutine typh_gather_z(in, dims, rank, out, typh_err)

    type(c_ptr), intent(in) :: in
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    type(c_ptr), intent(out) :: out
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    type(c_ptr) :: dims_c

    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_gather_z_c(in, dims_c, rank_c, out), kind=i32)

end subroutine typh_gather_z


subroutine typh_gather(datatype, in, dims, rank, out, typh_err)

    integer(kind(typh_datatype_null)), intent(in) :: datatype
    type(c_ptr), intent(in) :: in
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    type(c_ptr), intent(out) :: out
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    type(c_ptr) :: dims_c

    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_gather_c(datatype, in, dims_c, rank_c, out), kind=i32)

end subroutine typh_gather


subroutine typh_reduce_i(in, dims, rank, out, op, typh_err)

    type(c_ptr), intent(in) :: in
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    type(c_ptr), intent(out) :: out
    integer(kind(typh_op_sum)), intent(in) :: op
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    type(c_ptr) :: dims_c

    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_reduce_i_c(in, dims_c, rank_c, out, op), kind=i32)

end subroutine typh_reduce_i


subroutine typh_reduce_d(in, dims, rank, out, op, typh_err)

    type(c_ptr), intent(in) :: in
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    type(c_ptr), intent(out) :: out
    integer(kind(typh_op_sum)), intent(in) :: op
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    type(c_ptr) :: dims_c

    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_reduce_d_c(in, dims_c, rank_c, out, op), kind=i32)

end subroutine typh_reduce_d


subroutine typh_reduce_z(in, dims, rank, out, op, typh_err)

    type(c_ptr), intent(in) :: in
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    type(c_ptr), intent(out) :: out
    integer(kind(typh_op_sum)), intent(in) :: op
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    type(c_ptr) :: dims_c

    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_reduce_z_c(in, dims_c, rank_c, out, op), kind=i32)

end subroutine typh_reduce_z


subroutine typh_reduce(datatype, in, dims, rank, out, op, typh_err)

    integer(kind(typh_datatype_null)), intent(in) :: datatype
    type(c_ptr), intent(in) :: in
    integer(i32), intent(in) :: rank
    integer(i32), dimension(rank), target, intent(in) :: dims
    type(c_ptr), intent(out) :: out
    integer(kind(typh_op_sum)), intent(in) :: op
    integer(i32), intent(out) :: typh_err

    integer(c_int) :: rank_c
    type(c_ptr) :: dims_c

    rank_c = int(rank, kind=c_int)
    dims_c = c_loc(dims)

    typh_err = int(typh_reduce_c(datatype, in, dims_c, rank_c, out, op), &
        kind=i32)

end subroutine typh_reduce

end module typhon
