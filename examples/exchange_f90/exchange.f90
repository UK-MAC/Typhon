program exchange

use, intrinsic :: iso_c_binding, only: c_loc
use typhon

implicit none

    integer(i32), parameter :: ndat = 7

    integer(i32) :: typh_err

    integer(i32) :: i, j, ip, iel, ilay
    integer(i32) :: nproc, rank, sqrt_nproc, meshw, meshh, meshwl, meshhl, &
        xoff, yoff, nell, nndl, nelg, nndg, ndperel, nlay, global_el_num, &
        total_nel, total_nnd, partition_id, keyset_id, elquant_id, phase_id
    integer(i32), dimension(4) :: global_nd_num
    integer(i32), dimension(2) :: dims
    integer(i32), dimension(:), allocatable :: partition, layer_nel, &
        layer_nnd, layer_cum_nel, layer_cum_nnd
    integer(i32), dimension(:,:), allocatable :: conn_data
    real(f64), dimension(:,:), allocatable, target :: elquant
    integer(i32), pointer :: ellocglob(:), ndlocglob(:), elreg(:), &
        elmat(:), elnd(:,:), elowner(:,:), ndowner(:,:)
    logical(z32) :: is_master

    ! mesh dims
    meshw = 100
    meshh = 100


    !---------------------------------------------------------------------------
    ! INITIALISATION
    !---------------------------------------------------------------------------
    call typh_init(typh_err)
    call check_typh_err(typh_err, "typh_init")

    call typh_get_size(nproc, typh_err)
    call check_typh_err(typh_err, "typh_get_size")

    call typh_get_rank(rank, typh_err)
    call check_typh_err(typh_err, "typh_get_rank")

    call typh_is_master(is_master, typh_err)
    if (is_master) then
        print *, "nproc = ", nproc
    endif

    call typh_start_register(typh_err)
    call check_typh_err(typh_err, "typh_start_register")

    !
    !  O--                --> x=meshw
    !  |  x x x x . x x x x x
    !     x x x x . x x x x x
    !     x x x x . x x x x x
    !     x x x x . x x x x x
    !     x x x x . x x x x x
    !     . . . . . . . . . .
    !     x x x x . x x x x x
    !     x x x x . x x x x x
    !     x x x x . x x x x x
    !  |  x x x x . x x x x x
    !  v  x x x x . x x x x x
    ! y=meshh
    !
    ! Divide mesh into equal size squares. nproc must be a power of two and
    ! both mesh dims must divide evenly into sqrt(nproc)
    !
    sqrt_nproc = int(sqrt(real(nproc, kind=f64)), kind=i32)

    meshwl = meshw / sqrt_nproc
    meshhl = meshh / sqrt_nproc

    xoff = mod(rank, sqrt_nproc) * meshwl
    yoff = (rank / sqrt_nproc) * meshhl

    ! Print partitioning
    do ip=0,nproc
        if (rank == ip) then
            print *, "rank = ", rank, "(", xoff, ", ", yoff, ", ", &
                meshwl, ", ", meshhl, ")"
        endif
        call typh_barrier(typh_err)
        call check_typh_err(typh_err, "typh_barrier")
    end do

    ! Distribute mesh
    nell = meshwl * meshhl ! local # elements
    nndl = (meshwl+1) * (meshhl+1) ! local # nodes
    nelg = meshw * meshh ! global # elements
    nndg = (meshw+1) * (meshh+1) ! global # nodes
    ndperel = 4 ! square elements have 4 corners
    nlay = 1 ! 1 ghost layer

    if (is_master) then
        print *, "nelg = ", nelg
        print *, "nndg = ", nndg
    endif

    ! Set local element colour (each proc just starts with the elements it will
    ! ultimately control)
    allocate(partition(nell))
    partition(:) = rank

    ! Init connectivity data
    !
    !  - global element num
    !  - element region index
    !  - element material index
    !  - element global node num 1..4
    !
    allocate(conn_data(ndat, nell))
    iel = 1
    do i=1,meshwl
        do j=1,meshhl
            global_el_num = get_el_num(meshw, meshh, xoff+i-1, yoff+j-1)

            ! number nodes counter clockwise
            !
            ! 4     3
            !  x---x
            !  |   |
            !  x---x
            ! 1     2
            !
            global_nd_num(1) = get_nd_num(meshw, meshh, xoff+i-1, yoff+j  )
            global_nd_num(2) = get_nd_num(meshw, meshh, xoff+i  , yoff+j  )
            global_nd_num(3) = get_nd_num(meshw, meshh, xoff+i  , yoff+j-1)
            global_nd_num(4) = get_nd_num(meshw, meshh, xoff+i-1, yoff+j-1)

            conn_data(1, iel) = global_el_num
            conn_data(2, iel) = 0 ! don't care
            conn_data(3, iel) = 0 ! don't care
            conn_data(4:7, iel) = global_nd_num(:)
            iel = iel+1
        end do
    end do

    ! Call typh_distribute_mesh
    allocate(layer_nel(nlay+1))
    allocate(layer_nnd(nlay+1))
    total_nel = 0
    total_nnd = 0

    call typh_distribute_mesh(nell, nndl, nelg, nndg, ndperel, nlay, conn_data, &
        partition, layer_nel, layer_nnd, total_nel, total_nnd, ellocglob, &
        ndlocglob, elreg, elmat, elnd, elowner, ndowner, typh_err)
    call check_typh_err(typh_err, "typh_distribute_mesh")

    call typh_free(c_loc(elmat), typh_err)
    call typh_free(c_loc(elreg), typh_err)
    deallocate(conn_data)
    deallocate(partition)

    ! Init Typhon decomposition info
    ! ... Calculate cumulative layer sizes
    allocate(layer_cum_nel(nlay+1))
    allocate(layer_cum_nnd(nlay+1))
    layer_cum_nel(1) = layer_nel(1)
    layer_cum_nnd(1) = layer_nnd(1)
    do ilay=2,nlay+1
        layer_cum_nel(ilay) = layer_cum_nel(ilay-1) + layer_nel(ilay)
        layer_cum_nnd(ilay) = layer_cum_nnd(ilay-1) + layer_nnd(ilay)
    end do

    call typh_set_partition_info(partition_id, total_nel, total_nnd, &
        typh_shape_quad, nlay, &
        layer_cum_nel, layer_cum_nnd, elowner, ndowner, ellocglob, &
        ndlocglob, elnd, "partitioning", typh_err)
    call check_typh_err(typh_err, "typh_set_partition_info")

    deallocate(layer_cum_nnd)
    deallocate(layer_cum_nel)

    ! Create key set
    call typh_create_key_set(typh_keytype_cell, 1, nlay, partition_id, &
        keyset_id, typh_err)
    call check_typh_err(typh_err, "typh_create_key_set")

    ! Add a quant
    ! Set local values to our rank and ghosts to -1
    allocate(elquant(ndperel, total_nel))
    do iel=1,layer_nel(1)
        elquant(:,iel) = real(rank, kind=f64)
    end do

    do iel=layer_nel(1)+1,total_nel
        elquant(:,iel) = -1.0
    end do

    ! Register quant with Typhon
    ! XXX note we give typhon reversed dimensions, as it expects row-major
    !     layout and we are giving column-major
    dims(1) = typh_mesh_dim
    dims(2) = ndperel
    call typh_add_quant(elquant_id, "elquant", typh_ghosts_one, &
        typh_datatype_real, typh_centring_cell, .TRUE._z32, &
        typh_auxiliary_none, dims, 2, typh_err)
    call check_typh_err(typh_err, "typh_add_quant")

    if (elquant_id /= 0) then
        STOP 3
    end if

    dims(1) = total_nel
    call typh_set_quant_address(elquant_id, c_loc(elquant), dims, 2, typh_err)
    call check_typh_err(typh_err, "typh_set_quant_address")

    ! Add a phase
    call typh_add_phase(phase_id, "phase", typh_ghosts_one, .TRUE._z32, &
        keyset_id, -1, typh_err)
    call check_typh_err(typh_err, "typh_add_phase")

    if (phase_id /= 0) then
        STOP 2
    end if

    call typh_add_quant_to_phase(phase_id, elquant_id, -1, -1, -1, -1, &
        typh_err)
    call check_typh_err(typh_err, "typh_add_quant_to_phase")

    ! Initialisation complete
    call typh_finish_register(typh_err)
    call check_typh_err(typh_err, "typh_finish_register")


    !---------------------------------------------------------------------------
    ! EXCHANGE
    !---------------------------------------------------------------------------
    ! Before exchange, all ghosts are set to -1
    do iel=layer_nel(1)+1,total_nel
        do j=1,ndperel
            if (elquant(ndperel,iel) /= -1.0) then
                print *, "error at (", iel, ", ", j, ")"
                stop 1
            end if
        end do
    end do

    call typh_exchange(phase_id, typh_err)
    call check_typh_err(typh_err, "typh_exchange")

    ! After exchange, ghosts should contain their owner's rank
    do iel=layer_nel(1)+1,total_nel
        do j=1,ndperel
            if (elquant(ndperel,iel) /= real(elowner(1,iel), kind=f64)) then
                print *, "error at (", iel, ", ", j, ")"
                stop 1
            end if
        end do
    end do

    call typh_barrier(typh_err)
    call check_typh_err(typh_err, "typh_barrier")
    if (is_master) then
        print *, "exchange succeeded"
    end if


    !---------------------------------------------------------------------------
    ! CLEANUP
    !---------------------------------------------------------------------------
    deallocate(elquant)
    call typh_free(c_loc(ndowner), typh_err)
    call typh_free(c_loc(elowner), typh_err)
    call typh_free(c_loc(elnd), typh_err)
    call typh_free(c_loc(ndlocglob), typh_err)
    call typh_free(c_loc(ellocglob), typh_err)
    deallocate(layer_nnd)
    deallocate(layer_nel)

    call typh_kill(.TRUE._z32, typh_err)
    call check_typh_err(typh_err, "typh_kill")

contains

subroutine check_typh_err(typh_err, fname)
    integer(i32), intent(in) :: typh_err
    character(len=*), intent(in) :: fname

    if (typh_err /= typh_success) then
        print *, fname, "failed"
        stop 1
    endif
end subroutine check_typh_err

pure function get_el_num(meshw, meshh, x, y)
    integer(i32) :: get_el_num
    integer(i32), intent(in) :: meshw
    integer(i32), intent(in) :: meshh
    integer(i32), intent(in) :: x
    integer(i32), intent(in) :: y

    get_el_num = y * meshw + x
end function get_el_num

pure function get_nd_num(meshw, meshh, x, y)
    integer(i32) :: get_nd_num
    integer(i32), intent(in) :: meshw
    integer(i32), intent(in) :: meshh
    integer(i32), intent(in) :: x
    integer(i32), intent(in) :: y

    get_nd_num = y * (meshw+1) + x
end function get_nd_num

end program exchange
