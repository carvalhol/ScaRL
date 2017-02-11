program main_ScaRL

    use mpi
    use write_output
    use str_functions
    use fftw_ScaRL

    implicit none

    !INPUTS
    integer :: rank, nb_procs
    !LOCAL VARIABLES
    integer, dimension(3) :: Np=[5,10,15]
    integer :: comm_group
    integer, dimension(3) :: L=[20,30,40]
    integer, dimension(3) :: Np_ovlp = [3,4,5], topo_shape, topo_pos
    !INPUT VARIABLES

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, comm_group, rank, nb_procs)

    print *, "I'm processor number ", rank, " of ", nb_procs

    !Finding place in topology
    call decide_topo_shape(nb_procs, L, &
                            Np_ovlp, topo_shape)
    if(rank == 0) print *, "topo_shape =",topo_shape
    call get_topo_pos(rank, topo_shape, topo_pos)
    print *, "I'm processor number ", rank, "->topo_pos = ", topo_pos
    !Creating Fields
    call create_fields(Np(:), Np_ovlp, rank, nb_procs, &
                       topo_pos, topo_shape, comm_group)

    !Finalizing MPI
    call end_communication()



        !------------------------------------------------
        !------------------------------------------------
        !------------------------------------------------
        !------------------------------------------------
    contains

        !----------------------------------------------------------
        !----------------------------------------------------------
        !----------------------------------------------------------
        !----------------------------------------------------------
        subroutine init_communication(comm_local, comm, rank, nb_procs)
            implicit none
            !INPUT
            integer, intent(in) :: comm_local
            !OUTPUT
            integer, intent(out) :: comm, rank, nb_procs
            !LOCAL
            integer :: code

            call MPI_INIT(code)
            call MPI_COMM_RANK(comm_local, rank, code)
            call MPI_COMM_SIZE(comm_local, nb_procs, code)

            comm = comm_local

        end subroutine init_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine end_communication()
            implicit none
            !LOCAL
            integer :: code

            call MPI_FINALIZE(code)

        end subroutine end_communication

        !---------------------------------------------
        !---------------------------------------------
        !---------------------------------------------
        !---------------------------------------------
        subroutine create_fields(Np, Np_ovlp, rank, &
                                 nb_procs, topo_pos, topo_shape, &
                                 comm_group)
            implicit none
            !INPUT
            integer, dimension(3), intent(in) :: Np
            integer, dimension(3), intent(in) :: Np_ovlp
            integer, intent(in) :: rank, comm_group, nb_procs
            integer, dimension(3), intent(in) :: topo_pos, topo_shape
            double precision, dimension(3) :: xStep = [1d0, 1d0, 1d0] 
            double precision, dimension(3) :: xMinGlob = [-1d0,-2d0,-3d0]
             
            !LOCAL
            double precision, dimension(Np(1), Np(2), Np(3)) :: k_mtx
            character(len=1024) :: HDF5_name, XMF_name, res_folder
            integer, dimension(3) :: origin
            double precision, dimension(3) :: coord_0, coord_N
            double precision, dimension(3,nb_procs) :: coord_0_list
            character(len=1024), dimension(nb_procs) :: HDF5_list 
            integer :: temp_rank
            integer, dimension(3) :: temp_origin
            integer, dimension(3) :: temp_topo_pos
            logical :: oneFile=.true.
            integer :: partition_type = 1
            double precision, dimension(3) :: xRange
             
            res_folder = "."
            HDF5_name = str_cat("HDF5_proc_",trim(num2str(rank,4)),".h5")
            XMF_name  = str_cat("XMF_proc_",trim(num2str(rank,4)),".xmf")
            origin = (Np-Np_ovlp)*topo_pos            
            coord_0 = dble(origin)*xStep + xMinGlob
            coord_N = coord_0 + (dble(Np-1)*xStep)
            xRange = coord_N - coord_0

            !k_mtx(:,:,:) = dble(rank)
            !k_mtx(:,:,:) = 1d0

            call gen_std_gauss_FFT(k_mtx, Np, xRange)

            call apply_UnityPartition_mtx(Np, Np_ovlp,&
                                          partition_type, &
                                          topo_pos, topo_shape, &
                                          k_mtx)
            call add_overlap(k_mtx, Np, Np_ovlp, rank, &
                           nb_procs, topo_pos, topo_shape, &
                           comm_group)

            if(rank == 0) print*, "maxval(k_mtx) AFTER = ", maxval(k_mtx) 
            if(rank == 0) print*, "minval(k_mtx) AFTER = ", minval(k_mtx) 
            

            if(oneFile) then 
                call write_hdf5_multi_proc_3D(coord_0, &
                                          coord_N,     &
                                          k_mtx,       &
                                          "HDF5_global.h5",         &
                                          res_folder,           &
                                          rank, nb_procs,           &
                                          comm_group)
            else
                call write_hdf5_single_proc_3D(coord_0, &
                                           coord_N, &
                                           k_mtx, &
                                           HDF5_name, &
                                           XMF_name, &
                                           res_folder, &
                                           rank)
            end if

            if(rank==0) then
                do temp_rank = 0, nb_procs-1
                    call get_topo_pos(temp_rank,  &
                                      topo_shape, &
                                      temp_topo_pos)
                    HDF5_list(temp_rank+1) = &
                     str_cat("HDF5_proc_",trim(num2str(temp_rank,4)),".h5")
                    if(oneFile) HDF5_list(temp_rank+1) = "HDF5_global.h5" 
                    temp_origin = (Np-Np_ovlp)*temp_topo_pos            
                    coord_0_list(:, temp_rank+1) = &
                        dble(temp_origin)*xStep + xMinGlob
                end do
                call write_XMF_global("XMF_global.xmf", HDF5_list, &
                                 coord_0_list, xStep, &
                                 shape(k_mtx), nb_procs)
            end if

        end subroutine create_fields


        !-----------------------------------------------------------
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        subroutine get_topo_pos(rank, topo_shape, topo_pos)
        implicit none
        !INPUT
        integer, intent(in) :: rank
        integer, dimension(3), intent(in) :: topo_shape
        !OUTPUT
        integer, dimension(3), intent(out) :: topo_pos
        !LOCAL
        integer :: rest

        rest = rank
        topo_pos(3) = rest/product(topo_shape(1:2))
        rest = rest - topo_pos(3)*product(topo_shape(1:2))
        topo_pos(2) = rest/topo_shape(1)
        rest = rest - topo_pos(2)*topo_shape(1)
        topo_pos(1) = rest

        end subroutine get_topo_pos
        
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        subroutine get_topo_rank(topo_shape, topo_pos, rank)
        implicit none
        !INPUT
        integer, dimension(3), intent(in) :: topo_shape
        integer, dimension(3), intent(in) :: topo_pos
        !OUTPUT
        integer, intent(out) :: rank

        rank = product(topo_shape(1:2))*topo_pos(3) &
               + topo_shape(1)*topo_pos(2)          &
               + topo_pos(1)

        end subroutine get_topo_rank

        !-----------------------------------------------------------
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        subroutine decide_topo_shape(nb_procs, nPointsMtx, nPointsOvlp, topo_shape)

            implicit none
            !INPUT
            integer, intent(in) :: nb_procs
            integer, dimension(3), intent(in) :: nPointsMtx, nPointsOvlp

            !OUTPUT
            integer, dimension(3), intent(out) :: topo_shape
            !LOCAL
            integer, dimension(3) :: nPoints, nBlocks
            integer, dimension(3) :: nFieldsIdeal, nBlocksIdeal
            integer, dimension(3) :: nFieldsChosen
            integer, dimension(3) :: nPointsBase, ratio
            logical :: nFieldsOK
            integer, dimension(100) :: factors
            integer :: i, pos, np, np_total, np_start, np_end
            double precision :: vol_surf_factor, vol_surf_factor_temp
            integer :: pointsPerBlockMax = 10

            nPointsBase = nPointsMtx
            nBlocks = 1
            nFieldsOK = .false.

            do while (.not. nFieldsOK)
                nPoints = nPointsBase + (nBlocks - 1)*2*nPointsOvlp
                !print*, "nPoints = ", nPoints
                !print*, "nBlocks = ", nBlocks
                where(nPoints/nBlocks > pointsPerBlockMax) & 
                                                   nBlocks = nBlocks + 1
                ratio = ceiling(dble(nPoints)/dble(nBlocks))
                if(all(ratio < pointsPerBlockMax)) nFieldsOK = .true.
                if(product(nBlocks) > nb_procs) nFieldsOK = .true. 
            end do

            nBlocksIdeal = nBlocks


            np_total = nb_procs
            vol_surf_factor_temp = 0D0

            np_start = ceiling(0.9*np_total)
            np_end   = np_total

            do np = np_start, np_end
                nBlocks = nBlocksIdeal
                factors = 0
                nFieldsIdeal  = 1
                call find_factors(np, factors)

                do i = size(factors), 1, -1
                    if(factors(i) == 0) cycle
                    if(all(nBlocks == 1)) exit
                    pos = MAXLOC(nBlocks,1)
                    nFieldsIdeal(pos) = nFieldsIdeal(pos)*factors(i)
                    nBlocks(pos) = &
                         ceiling(dble(nBlocks(pos))/dble(factors(i)))
                end do

                vol_surf_factor_temp = &
                   dble(product(nFieldsIdeal))/&
                   dble(2*sum(nFieldsIdeal*CSHIFT(nFieldsIdeal, shift=1)))

                if(vol_surf_factor_temp >= vol_surf_factor) then
                    vol_surf_factor = vol_surf_factor_temp
                    nFieldsChosen   = nFieldsIdeal
                end if
            end do
            topo_shape = nFieldsChosen

        end subroutine decide_topo_shape

        !-----------------------------------------------------------
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        !-----------------------------------------------------------
        subroutine find_factors(n, d)
        implicit none 
        integer, intent(in) :: n
        integer, dimension(:), intent(out) :: d
        integer :: div, next, rest
        integer :: i

        i = 1
        div = 2; next = 3; rest = n

        do while ( rest /= 1 )
           do while ( mod(rest, div) == 0 )
              d(i) = div
              i = i + 1
              rest = rest / div
           end do
           div = next
           next = next + 2
        end do
        end subroutine find_factors

    !--------------------------------------------------------------
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    !--------------------------------------------------------------
    subroutine add_overlap(RF, Np, Np_ovlp, rank, &
                           nb_procs, topo_pos, topo_shape, &
                           comm_group)
        implicit none
        !INPUT
        integer, dimension(3), intent(in) :: Np
        integer, dimension(3), intent(in) :: Np_ovlp
        integer, intent(in) :: rank, comm_group, nb_procs
        integer, dimension(3), intent(in) :: topo_pos, topo_shape

        !INPUT OUTPUT
        double precision, dimension(:,:,:), intent(inout) :: RF

        !LOCAL
        double precision, dimension(Np(1), Np(2), Np(3)) :: RF_temp1, RF_temp2
        integer, dimension(27) :: neigh_rank, op_neigh_rank
        integer, dimension(3,27) :: neigh_shift
        integer, dimension(-1:1,-1:1,-1:1) :: dir_index
        integer, dimension(3) :: topo_pos_temp
        !LOCAL
        integer :: dir, op_dir
        integer, dimension(3) :: minP, maxP
        integer :: neighRank
        integer, dimension(3) :: dirShift
        integer(kind=8) :: nOvlpMax
        integer :: code
        integer :: double_size
        integer(kind=8) :: totalSize, overHead, overEst, bufferSize
        double precision, dimension(:), allocatable :: buffer
        integer :: request
        integer, dimension(MPI_STATUS_SIZE) :: statut
        integer :: tag
        logical :: snd, rcv
        integer :: BufDT_size
        integer :: ii, jj, kk, cnt
        integer :: rang_test=0

        !Defining neighbours
        cnt = 0
        do kk = -1, 1
            do jj = -1, 1
                do ii = -1, 1
                    cnt = cnt+1
                    neigh_shift(:,cnt) = [ii,jj,kk]
                    dir_index(ii,jj,kk) = cnt
                    
                    !Direct neighbour
                    topo_pos_temp = topo_pos+[ii,jj,kk]
                    if(any(topo_pos_temp < 0) .or. any(topo_pos_temp > (topo_shape -1))) then
                        neigh_rank(cnt) = rank 
                    else
                       call  get_topo_rank(topo_shape, &
                                           topo_pos_temp, &
                                            neigh_rank(cnt)) 
                    end if
                    
                    !Oposite neighbour
                    topo_pos_temp = topo_pos-[ii,jj,kk]
                    if(any(topo_pos_temp < 0) .or. any(topo_pos_temp > (topo_shape -1))) then
                        op_neigh_rank(cnt) = rank 
                    else
                       call  get_topo_rank(topo_shape, &
                                           topo_pos_temp, &
                                           op_neigh_rank(cnt)) 
                    end if
                    
                   ! if(neigh_rank(cnt) /= rank .or. op_neigh_rank(cnt) /= rank)then 
                   ! if(rang_test == rank) print*, "rk ", rank, "neigh shift", neigh_shift(:,cnt), "---------------"
                   ! if(rang_test == rank .and. neigh_rank(cnt) /= rank) &
                   !       print*, "neigh = ", neigh_rank(cnt)
                   ! if(rang_test == rank .and. op_neigh_rank(cnt) /= rank) &
                   !       print*, "op_ng = ", op_neigh_rank(cnt) 
                   ! if(rang_test == rank) print*, "  " 
                   ! end if

                end do
            end do
        end do

        !Set Max Number of Points on Overlap (define buffer size)
        nOvlpMax = 0
        do dir = 1, size(neigh_rank)
            if(neigh_rank(dir) == rank) cycle 
            minP(:) = 1
            maxP(:) = Np
            where(neigh_shift(:,dir) == 1 ) minP = Np - Np_ovlp + 1
            where(neigh_shift(:,dir) == -1) maxP = Np_ovlp
            nOvlpMax = nOvlpMax + product(maxP - minP +1)
        end do

        !Buffer allocation
        overEst = 2
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,double_size,code)
        overHead = int(1+(MPI_BSEND_OVERHEAD)/double_size)
        bufferSize = overEst*(nOvlpMax+overHead)
        allocate(buffer(bufferSize))
        call MPI_BUFFER_ATTACH(buffer, int(double_size*bufferSize),code)

        !Communications
        RF_temp1 = 0

        do dir = 1, size(neigh_rank)

            call MPI_BARRIER(comm_group, code) !This is necessary so we can reduce the size of the buffer

            !SENDING---------------------------------------------------------------

            snd = .true.
            dirShift   = neigh_shift(:, dir)
            neighRank  = neigh_rank(dir)
            minP(:) = 1
            maxP(:) = Np
            where(neigh_shift(:,dir) == 1 ) minP = Np - Np_ovlp + 1
            where(neigh_shift(:,dir) == -1) maxP = Np_ovlp
            totalSize = (product(maxP - minP + 1))

            if(neigh_rank(dir) == rank) snd = .false. !Check if this direction exists

            if(snd) then
                !call wLog(" SENDING ==============")
                
                tag = neighRank

                call MPI_IBSEND (RF(minP(1):maxP(1), &
                                           minP(2):maxP(2), &
                                           minP(3):maxP(3)), &
                        int(totalSize), MPI_DOUBLE_PRECISION, &
                        neighRank, tag, comm_group, request, code)
            else

            end if

            !RECEIVING---------------------------------------------------------------
            
            rcv = .true.
            !dirShift   = -dirShift
            neighRank  = op_neigh_rank(dir)
            op_dir = dir_index(-dirShift(1), -dirShift(2), -dirShift(3))
            minP(:) = 1
            maxP(:) = Np
            where(neigh_shift(:,op_dir) == 1 ) minP = Np - Np_ovlp + 1
            where(neigh_shift(:,op_dir) == -1) maxP = Np_ovlp
            totalSize = (product(maxP - minP + 1))

            if(neigh_rank(op_dir) == rank) rcv = .false. !Check if this direction exists

            if(rcv) then
                !if(rank == rang_test) print*, " RECEIVING =============="

                tag  = rank
                RF_temp2 = 0.0D0

                call MPI_RECV (RF_temp2(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3)), &
                    int(totalSize), MPI_DOUBLE_PRECISION, &
                    neighRank, tag, comm_group, statut, code)
                    
                RF_temp1(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3)) = &
                    RF_temp1(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))   &
                    + RF_temp2(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))
                
               ! if(rank == rang_test) print*, "  RECVER: rank ", rank
               ! if(rank == rang_test) print*, "  FROM rank ", neighRank
               ! if(rank == rang_test) print*, "  minP = ", minP
               ! if(rank == rang_test) print*, "  maxP = ", maxP
               ! if(rank == rang_test) print*, "  op_dir_shft = ", -dirShift
               ! if(rank == rang_test) print*, "  totalSize = ", totalSize
               ! if(rank == rang_test) print*, "  tag = ", tag
               ! if(rank == rang_test) print*, "  CONTENT = ", RF_temp2(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))
            else
            
            end if
        end do

        RF = RF + RF_temp1

        BufDT_size = int(double_size*bufferSize)
        call MPI_BUFFER_DETACH (buffer,BufDT_size,code)
        
        if(allocated(buffer)) deallocate(buffer)

    end subroutine add_overlap

    !---------------------------------------------------------------
    !---------------------------------------------------------------
    !---------------------------------------------------------------
    !---------------------------------------------------------------
    subroutine apply_UnityPartition_mtx(Np, Np_ovlp,&
                                        partitionType, &
                                        topo_pos, topo_shape, &
                                        randField_3D)
    
        implicit none
    
        !INPUT
        integer, dimension(3), intent(in) ::  Np, Np_ovlp
        integer, intent(in) :: partitionType
        integer, dimension(3) :: topo_pos, topo_shape
        double precision, dimension(:,:,:) :: randField_3D
        !LOCAL
        integer :: ii, jj, kk
        double precision, dimension(:), allocatable :: pattern
        integer, dimension(3) :: U_Lim, D_Lim
        logical :: considerNeighbour
        double precision, dimension(Np(1), Np(2), Np(3)) :: unityPartition
        double precision, parameter :: PI = 3.1415926535898d0;
        integer :: shift

        unityPartition = 1.0D0
    
        !print*, "apply_UnityPartition_mtx = ", pattern
        !print*, "topo_pos = ", topo_pos
        !print*, "topo_shape = ", topo_shape
        
        do jj = 1, 2
    
            if(jj == 1) then
                D_Lim = 1
                U_Lim = Np_ovlp
                shift = -1
            else
                D_Lim = Np - Np_ovlp + 1
                U_Lim = Np
                shift = 1
            end if
    
            do ii = 1, 3
    
                if(topo_pos(ii) + shift < 0 &
                   .or. topo_pos(ii) + shift > (topo_shape(ii) - 1)) cycle 
    
                if(allocated(pattern)) deallocate(pattern)
                allocate(pattern(Np_ovlp(ii)))
                pattern = 1d0
                if(Np_ovlp(ii) > 1) then
                   pattern = [(dble(kk), kk=0, Np_ovlp(ii)-1)]/dble((Np_ovlp(ii)-1))
                   if(jj == 1) pattern = pattern(size(pattern):1:-1) !Pattern inversion
                   if(partitionType == 1) then
                       pattern = ((1.0D0 + cos(PI*(pattern)))&
                                   / 2.0D0)
                   end if
                end if
                !print*, "pattern = ", pattern
                !print*, "size(pattern) = ", size(pattern)
                !print*, "D_Lim = ", D_Lim
                !print*, "U_Lim = ", U_Lim
    
                if(ii==1) then
                    do kk = D_Lim(ii), U_Lim(ii)
                        unityPartition(kk,:,:) = unityPartition(kk,:,:) & 
                                                 * pattern(kk-D_Lim(ii)+1)
                    end do
                else if(ii==2) then
                    do kk = D_Lim(ii), U_Lim(ii)
                        unityPartition(:,kk,:) = unityPartition(:,kk,:) &
                                                 * pattern(kk-D_Lim(ii)+1)
                    end do
                else if(ii==3) then
                    do kk = D_Lim(ii), U_Lim(ii)
                        unityPartition(:,:,kk) = unityPartition(:,:,kk) &
                                                 * pattern(kk-D_Lim(ii)+1)
                    end do
                end if
    
            end do
        end do
    
        randField_3D = randField_3D * sqrt(unityPartition)
    
        if(allocated(pattern)) deallocate(pattern)
    
    end subroutine apply_UnityPartition_Mtx
end program main_ScaRL



!module fftw_ScaRL
!    use, intrinsic :: iso_c_binding
!    !include 'fftw3.f'
!    include 'fftw3-mpi.f03'
!    !Attention there is a line "include 'fftw3-mpi.f03'" just before this one
!
!contains
!
!!-------------------------------------------------------------
!!-------------------------------------------------------------
!!-------------------------------------------------------------
!!-------------------------------------------------------------
!
!
!subroutine gen_std_gauss_FFT(data_real_3D, Np, xRange)
!
!    implicit none
!    !INPUT
!    integer, dimension(3), intent(in) :: Np
!    double precision, dimension(3), intent(in) :: xRange
!    integer :: corrMod=1
!    double precision, dimension(3) :: corrL=[1.0,1.0,1.0]
!    !INPUT OUTPUT
!    double precision, dimension(:,:,:), intent(inout) :: data_real_3D
!    !LOCAL
!    integer :: L, M, N
!    !integer(C_INTPTR_T) :: L, M, N
!    !integer(kind=8) :: plan
!    type(C_PTR) :: plan
!    double precision, dimension(Np(1),Np(2),Np(3)) :: phiK, gammaK, Sk_mtx
!    double precision :: ampMult
!    double precision :: k_max = 7.355d0
!    double precision, dimension(3) :: delta_k
!    integer :: ii, jj, kk
!    double precision, parameter :: PI = 3.1415926535898d0;
!    integer, parameter :: cm_GAUSSIAN = 1
!    !real(C_DOUBLE), pointer :: data_C_3D(:,:,:)
!
!    L = Np(1)
!    M = Np(2)
!    N = Np(3)
!
!    !Build Sk matrix
!    delta_k(:) = k_max/Np(:)
!    select case(corrMod)
!
!    case(cm_GAUSSIAN)
!    do ii = 1, Np(1)
!        do jj = 1, Np(2)
!            do kk = 1, Np(3)
!                Sk_mtx(ii,jj,kk) = product(corrL) &
!                       *exp(-sum((dble([ii,jj,kk]-1)*delta_k(:))**2.0D0 &
!                       *(corrL(:)**2.0D0)/(4.0d0*PI)))
!            end do
!        end do
!    end do
!
!    end select
!    
!    call random_number(gammaK(:,:,:))
!    call random_number(phiK(:,:,:))
!    gammaK  = gammaK -0.5
!    Sk_mtx  = gammak*sqrt(Sk_mtx)*cos(2.0D0*PI*phik);
!    
!    plan = fftw_plan_r2r_3d(N,M,L, Sk_mtx, data_real_3d, &
!           FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
!    call fftw_execute_r2r(plan, Sk_mtx, data_real_3D)
!    call fftw_destroy_plan(plan)
!
!end subroutine gen_Std_Gauss_FFT
!
!
!end module fftw_ScaRL
