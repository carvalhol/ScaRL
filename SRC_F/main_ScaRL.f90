program main_ScaRL

    use mpi
    use write_output
    use read_input
    use str_functions
    use fftw_ScaRL
    use constants_ScaRL

    implicit none

    !INPUTS SYSTEM
    integer :: rank, nb_procs
    !INPUTS FILE
    integer :: nSamples = 3
    character (len=1024) :: output_folder = "SAMPLES"
    double precision, dimension (:,:), allocatable :: xMinGlob
    double precision, dimension (:,:), allocatable :: xMaxGlob
    double precision, dimension (:,:), allocatable :: corrL
    integer, dimension(:), allocatable :: corrMod, margiFirst
    integer, dimension(:), allocatable :: seedBase
    character (len=1024), dimension(:), allocatable :: output_name
    double precision, dimension(:), allocatable :: avg, std_dev
    double precision, dimension (:,:), allocatable :: overlap
    integer, dimension(:,:), allocatable :: pointsPerCorrL
    
    !INPUTS HARD
    integer(kind = 8) :: pointsPerBlockIdeal = 100*100*100
       
    !LOCAL VARIABLES
    integer, dimension(3) :: Np, Np_temp
    integer :: comm_group
    integer, dimension(3) :: L
    integer, dimension(3) :: Np_ovlp, topo_shape, topo_pos
    double precision, dimension(3) :: xStep
    integer :: color, nb_procs_tmp, rank_tmp, code
    character (len=1024) :: command
    integer :: s

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, comm_group, rank, nb_procs)

    call read_input_ScaRL("./input_ScaRL.txt", rank, &
                          nSamples, output_folder, &
                          xMinGlob, xMaxGlob, &
                          corrL, corrMod, margiFirst, &
                          seedBase, &
                          output_name, &
                          avg, std_dev, overlap, &
                          pointsPerCorrL, comm_group)
    
    
    !Creating Output Folder 
    if(rank == 0) then
        command = 'mkdir -pv '// trim(adjustL(output_folder)) 
        call system(command)
    end if
    call MPI_BARRIER (comm_group ,code)

    do s = 1, nSamples !----------------------------
    if(rank == 0) print*, "==================================================== "
    if(rank == 0) print*, "Generating SAMPLE ",s 
    if(rank == 0) print*, " "
    if(rank == 0) print*, " "
    if(rank == 0) print*, " "
    !output_name(s) = str_cat("sample_", num2str(s))
    !avg(s) = dble(s)
    !CV(s) = dble(s)
    !seedBase(s) = s
    !xMinGlob(:,s) = [0d0, 0d0, 0d0]
    !xMaxGlob(:,s) = [10d0, 10d0, 10d0]
    !corrL(:,s)    = [1d0, 2d0, 3d0]
    !corrMod(s) = cm_GAUSSIAN
    !margiFirst(s) = fom_LOGNORMAL
    !overlap(:,s) = [5d0,5d0,5d0]
    !pointsPerCorrL(:,s) = [5,5,5]
     
    call verify_inputs(xMinGlob(:,s), xMaxGlob(:,s), corrL(:,s), & 
                       overlap(:,s), avg(s), std_dev(s), corrMod(s), margiFirst(s),&
                       pointsPerCorrL(:,s), pointsPerBlockIdeal)

    !Defining L, Np_ovlp and xStep (Np stands for Number of points)
    L       = 1+ceiling((xMaxGlob(:,s)-xMinGlob(:,s))/corrL(:,s))*(pointsPerCorrL(:,s)-1)
    Np_ovlp = ceiling(overlap(:,s)*dble(pointsPerCorrL(:,s)))
    xStep   = corrL(:,s)/(dble(pointsPerCorrL(:,s)-1))  

    !Finding place in topology and defining the number of points per processor (Np)
    call decide_topo_shape(nb_procs, L, Np_ovlp, &
                           pointsPerBlockIdeal, & 
                           topo_shape)
    
    !Changing to keep only procs used in this topology
    nb_procs_tmp = nb_procs
    nb_procs = product(topo_shape)
    color = 0
    if(rank >= nb_procs) color = 1
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, rank, comm_group, code)
    if(nb_procs_tmp /= nb_procs) then
        if(rank == 0) print*, "WARNING!! Ignoring ", nb_procs_tmp-nb_procs ,"procs"
        if(rank == 0) print*, "          Consider using ", nb_procs ,"procs next time"
    end if

    if(color == 0) then
        call MPI_COMM_SIZE(comm_group, nb_procs_tmp, code)        
        call MPI_COMM_RANK(comm_group, rank_tmp, code)
        if(nb_procs_tmp /= nb_procs .or. rank_tmp /= rank) &
            stop("ERROR!! Communicator ill-defined")

        if(std_dev(s) == 0d0) then
            if(rank == 0) then
                print*, "constant field"
                call write_constant_field(xMinGlob(:,s), &
                                                     xMaxGlob(:,s), &
                                                     avg(s),        &
                                                     output_name(s), &
                                                     output_folder)
            end if
            cycle
        end if

        Np = ceiling(dble(L + ((topo_shape-1)*Np_ovlp))/dble(topo_shape))

        !Condition to not overlap more than 8 overlap areas
        Np_temp = Np
        where(topo_shape == 2 .and. Np < Np_ovlp)   Np = Np_ovlp
        where(topo_shape >= 3 .and. Np < 2*Np_ovlp) Np = 2*Np_ovlp
        if(any(Np_temp /= Np)) then
            if(rank == 0) print *, "WARNING!! Overlap is too big for this domain"
            if(rank == 0) print *, "                  the domain will be expanded"
            if(rank == 0) print *, "        (Consider using a smaller &&
                                    overlap/domain size (ratio) or fewer procs)"
            if(rank == 0) print *, "BEFORE: Np = ",  Np_temp, "Np_ovlp = ", Np_ovlp
            if(rank == 0) print *, "AFTER: Np = ",  Np, "Np_ovlp = ", Np_ovlp
        end if

        L = Np*topo_shape - ((topo_shape-1)*Np_ovlp)
        if(rank == 0) print *, "nb_procs   =", nb_procs
        if(rank == 0) print *, "topo_shape =", topo_shape
        if(rank == 0) print *, "L          =", L
        if(rank == 0) print *, "Np         =", Np
        if(rank == 0) print *, "Np_ovlp    =", Np_ovlp
        if(rank == 0) print *, "product(Np)=", product(Np)
        if(rank == 0) print *, "xStep      =", xStep

        call get_topo_pos(rank, topo_shape, topo_pos)
        

        !Creating Fields
        call create_fields(Np, Np_ovlp, L, pointsPerCorrL(:,s), rank, &
                           nb_procs, topo_pos, topo_shape, &
                           xStep, xMinGlob(:,s), &
                           corrL(:,s), corrMod(s), &
                           seedBase(s), & 
                           margiFirst(s), avg(s), std_dev(s), &
                           comm_group, &
                           output_name(s), &
                           output_folder)
    end if

    if(rank == 0) print*, " "
    if(rank == 0) print*, " "
    if(rank == 0) print*, "END Generating SAMPLE ",s 
    if(rank == 0) print*, "==================================================== "
    if(rank == 0) print*, " "
    if(rank == 0) print*, " "
    end do !---------------------------------


    if(allocated(output_name)) deallocate(output_name)

    if(rank == 0) print*, "EXIT CODE: OK " 
    
    !Finalizing MPI
    call MPI_COMM_FREE (comm_group, code)
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

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine verify_inputs(xMinGlob, xMaxGlob, corrL, &
                                 overlap, avg, std_dev, corrMod, margiFirst, & 
                                 pointsPerCorrL, pointsPerBlockIdeal)
       
        implicit none 
        double precision, dimension (3), intent(in) :: xMinGlob
        double precision, dimension (3), intent(in) :: xMaxGlob
        double precision, dimension (3), intent(in) :: corrL
        double precision, dimension (3), intent(in) :: overlap
        double precision, intent(in) :: avg
        double precision, intent(in) :: std_dev
        integer, intent(in) :: corrMod, margiFirst
        integer, dimension(3), intent(in) :: pointsPerCorrL
        integer(kind=8), intent(in) :: pointsPerBlockIdeal

        if(any(xMinGlob >= xMaxGlob)) then
            print*, "xMinGlob = ", xMinGlob
            print*, "xMaxGlob = ", xMaxGlob
            stop("ERROR: xMinGlob >= xMaxGlob")
        end if 

        if(any(corrL < 0d0)) then
            print*, "corrL = ", corrL
            stop("ERROR: corrL < 0d0")
        end if 
        
        if(any(overlap < 0d0)) then
            print*, "overlap = ", corrL
            print*, "WARNING: overlap < 0d0, the mesh will be discontinuous"
        end if 

        if((corrMod /= cm_GAUSSIAN)) then
            print*, "corrMod = ", corrMod
            stop("ERROR: Correlation model not implemented")
        end if
         
        if((margiFirst /= fom_GAUSSIAN) .and. &
            (margiFirst /= fom_LOGNORMAL)) then
            print*, "margiFirst = ", margiFirst
            stop("ERROR: First-order marginal not implemented")
        end if
        
        if((avg <= 0d0) .and. (margiFirst == fom_LOGNORMAL)) then
            print*, "       avg = ", avg
            print*, "margiFirst = ", margiFirst
            stop("ERROR: When using lognormal first-order marginal density the average should be greater than 0")
        end if 
        
        if((std_dev < 0d0)) then
            print*, " std_dev = ", std_dev
            stop("ERROR: Standard deviation should not be smaller than 0d0")
        end if
         
        if(any(pointsPerCorrL < 2)) then
            print*, "pointsPerCorrL = ", pointsPerCorrL
            stop("ERROR: PointsPerCorrL < 2")
        end if 
        
        if(pointsPerBlockIdeal < 8) then
            print*, "pointsPerBlockIdeal = ", pointsPerBlockIdeal
            stop("ERROR: pointsPerBlockIdeal < 8")
        end if 
        
        end subroutine verify_inputs
        !---------------------------------------------
        !---------------------------------------------
        !---------------------------------------------
        !---------------------------------------------
        subroutine create_fields(Np, Np_ovlp, L, pointsPerCorrL,  rank, &
                                 nb_procs, topo_pos, topo_shape, &
                                 xStep, xMinGlob, &
                                 corrL, corrMod, &
                                 seedBase, &
                                 margiFirst, avg, std_dev, &
                                 comm_group, &
                                 output_name, res_folder)
            implicit none
            !INPUT
            integer, dimension(3), intent(in) :: Np
            integer, dimension(3), intent(in) :: Np_ovlp
            integer, dimension(3), intent(in) :: L
            integer, dimension(3), intent(in) :: pointsPerCorrL
            integer, intent(in) :: rank, comm_group, nb_procs
            integer, dimension(3), intent(in) :: topo_pos, topo_shape
            double precision, dimension(3), intent(in) :: xStep
            double precision, dimension(3), intent(in) :: xMinGlob
            double precision, dimension(3), intent(in) :: corrL
            integer, intent(in) :: corrMod, margiFirst
            integer, intent(in) :: seedBase
            double precision, intent(in) :: avg, std_dev
            character(len=*), intent(in) :: output_name, res_folder
             
            !LOCAL
            double precision, dimension(Np(1), Np(2), Np(3)) :: k_mtx
            character(len=1024) :: HDF5_name, XMF_name
            double precision, dimension(3) :: coord_0, coord_N
            integer         , dimension(3) :: pos_0, pos_N
            double precision, dimension(3) :: xMaxGlob
            double precision, dimension(3,nb_procs) :: coord_0_list
            character(len=1024), dimension(nb_procs) :: HDF5_list 
            integer :: temp_rank
            integer, dimension(3) :: temp_origin
            integer, dimension(3) :: temp_topo_pos
            !logical :: oneFile=.false.
            logical :: oneFile=.true., oneDataSet=.true.
            integer :: partition_type = 1
            integer :: seedStart
            double precision, dimension(3) :: xRange, overlap
            integer :: code
            
            pos_0 = (Np-Np_ovlp)*topo_pos + 1
            pos_N = pos_0 + Np-1            
            coord_0 = dble(pos_0-1)*xStep + xMinGlob
            coord_N = dble(pos_N-1)*xStep + xMinGlob
            xRange = coord_N - coord_0
            overlap = dble(Np_ovlp - 1)*xStep 

            !Discover xMaxGlob
            call get_topo_pos(nb_procs-1,  &
                              topo_shape, &
                              temp_topo_pos)
            temp_origin = (Np-Np_ovlp)*temp_topo_pos            
            xMaxGlob = dble(temp_origin)*xStep + xMinGlob + xRange
            

            !k_mtx(:,:,:) = dble(rank)
            k_mtx(:,:,:) = 1d0
            if(rank == 0) print *, "Np                  = ", Np
            if(rank == 0) print *, "Np_ovlp             = ", Np_ovlp
            if(rank == 0) print *, "topo_pos            = ", topo_pos
            if(rank == 0) print *, "pos_0               = ", pos_0
            if(rank == 0) print *, "pos_N               = ", pos_N
            if(rank == 0) print *, "coord_0             = ", coord_0
            if(rank == 0) print *, "coord_N             = ", coord_N
            if(rank == 0) print *, "xMinGlob            = ", xMinGlob
            if(rank == 0) print *, "xMaxGlob            = ", xMaxGlob
            if(rank == 0) print *, "xRange (processor)  = ", xRange
            if(rank == 0) print *, "xRange (total)      = ", xMaxGlob - xMinGlob

            seedStart = seedBase + 1
            if(seedBase < 0) then
                !Stochastic Seed
                if(rank == 0) call calculate_random_nb(seedStart)
                call MPI_BCAST (seedStart, 1, MPI_INTEGER, 0, comm_group, code)
                if(seedStart == 0) seedStart = seedStart+1
                if(seedStart < 0) seedStart = abs(seedStart)  
            end if
            seedStart = seedStart + rank
            !print*, "Proc = ", rank, "seedStart = ", seedStart
            
            if(rank == 0) print*, "gen_std_gauss_FFT " 
            !call gen_std_gauss_FFT(k_mtx, Np, &
            !                       xRange, corrL, corrMod, seedStart)
            call gen_std_gauss_Shino_FFT(k_mtx, Np, &
                                   xRange, corrL, &
                                   pointsPerCorrL, corrMod, &
                                   seedStart, rank)

            !if(rank == 0) print*, "apply_UnityPartition " 
            !call apply_UnityPartition_mtx(Np, Np_ovlp,&
            !                              partition_type, &
            !                              topo_pos, topo_shape, &
            !                              k_mtx)
            !
            !if(rank == 0) print*, "add_overlap " 
            !call add_overlap(k_mtx, Np, Np_ovlp, rank, &
            !               nb_procs, topo_pos, topo_shape, &
            !               comm_group)
            
            !if(rank == 0) print*, "normalize_field " 
            !call normalize_field(topo_pos, topo_shape, Np, Np_ovlp, &
            !                     rank, comm_group, k_mtx)
            
            !if(rank == 0) print*, "multivariateTransformation " 
            !call multiVariateTransformation(avg, std_dev, margiFirst, &
            !                                k_mtx)
            
            if(rank == 0) print*, "maxval(k_mtx) AFTER = ", maxval(k_mtx) 
            if(rank == 0) print*, "minval(k_mtx) AFTER = ", minval(k_mtx)

            if(rank == 0) print*, "write HDF5 " 
            if(oneFile) then 
                if(oneDataSet) then
                    call write_hdf5_multi_proc_3D_1ds(coord_0, &
                                     coord_N,     &
                                     xMinGlob, xMaxGlob, &
                                     pos_0, pos_N, &
                                     L, Np, Np_ovlp, &
                                     k_mtx,       &
                                     str_cat(output_name,".h5"),  &
                                     str_cat(output_name,".xmf"),  &
                                     res_folder,           &
                                     rank, nb_procs,           &
                                     comm_group)
                    if(rank == 0) then
                        call write_HDF5_attributes(str_cat(res_folder,"/",output_name,".h5"), &
                                     nb_procs, 3, 1, FFT, seedBase, &
                                     corrMod, margiFirst, &
                                     topo_shape, &
                                     xMinGlob, xMaxGlob, xStep, corrL, overlap, &
                                     .false.)
                    end if
                else 
                    call write_hdf5_multi_proc_3D(coord_0, &
                                     coord_N,     &
                                     xMinGlob, xMaxGlob, &
                                     pos_0, pos_N, &
                                     k_mtx,       &
                                     str_cat(output_name,".h5"),  &
                                     res_folder,           &
                                     rank, nb_procs,           &
                                     comm_group)
                end if
            else
                HDF5_name = str_cat(output_name,"_proc_",trim(num2str(rank,5)),".h5")
                XMF_name  = str_cat(output_name,"_proc_",trim(num2str(rank,5)),".xmf")
                call write_hdf5_single_proc_3D(coord_0, &
                                           coord_N, &
                                           xMinGlob, xMaxGlob, &
                                           pos_0, pos_N, &
                                           k_mtx, &
                                           HDF5_name, &
                                           XMF_name, &
                                           res_folder, &
                                           rank)
            end if

            if(rank==0 .and. (.not.(oneFile .and.oneDataSet))) then
                if(rank == 0) print*, "WRITING GLOBAL XMF" 
                do temp_rank = 0, nb_procs-1
                    call get_topo_pos(temp_rank,  &
                                      topo_shape, &
                                      temp_topo_pos)
                    HDF5_list(temp_rank+1) = &
                     str_cat(output_name,"_proc_",trim(num2str(temp_rank,5)),".h5")
                    if(oneFile)  HDF5_list(temp_rank+1) = &
                                  str_cat(output_name,".h5") 
                    temp_origin = (Np-Np_ovlp)*temp_topo_pos            
                    coord_0_list(:, temp_rank+1) = &
                        dble(temp_origin)*xStep + xMinGlob
                end do
                call write_XMF_global(str_cat(res_folder,"/",output_name,".xmf"), &
                                 HDF5_list, &
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
        subroutine decide_topo_shape(nb_procs, L, nPointsOvlp, &
                                     pointsPerBlockIdeal, & 
                                     topo_shape)

            implicit none
            !INPUT
            integer, intent(in) :: nb_procs
            integer, dimension(3), intent(in) :: L, nPointsOvlp
            integer(kind = 8), intent(in) :: pointsPerBlockIdeal

            !OUTPUT
            integer, dimension(3), intent(out) :: topo_shape
            !LOCAL
            integer, dimension(3) :: nPoints, nBlocks
            integer, dimension(3) :: nFieldsIdeal, nBlocksIdeal
            integer, dimension(3) :: nFieldsChosen
            integer, dimension(3) :: nPointsBase, nPointsPerProc, ratio
            integer(kind = 8) :: nTotalPointsProc
            logical :: nFieldsOK
            integer, dimension(100) :: factors
            integer :: i, pos, np, np_total, np_start, np_end
            double precision :: vol_surf_factor, vol_surf_factor_temp

            nPointsBase = L
            nBlocks = 1
            nFieldsOK = .false.
            do while (.not. nFieldsOK)
            
                nPoints = nPointsBase + (nBlocks - 1)*nPointsOvlp
                nPointsPerProc = ceiling(dble(nPoints)/dble(nBlocks))
                nTotalPointsProc = product(int(nPointsPerProc,8))
                if(nTotalPointsProc <= pointsPerBlockIdeal) then
                    nFieldsOK = .true.
                else
                    pos = maxloc(nPointsPerProc, 1)
                    nBlocks(pos) = nBlocks(pos) + 1
                    if(product(nBlocks) > nb_procs) nFieldsOK = .true. 
                end if

                !print*, "nPoints = ", nPoints
                !print*, "nPointsPerProc = ", nPointsPerProc
                !print*, "nBlocks = ", nBlocks
            
            end do

            nBlocksIdeal = nBlocks


            np_total = nb_procs
            vol_surf_factor_temp = 0D0

            np_start = ceiling(0.75*np_total)
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
        subroutine normalize_field(topo_pos, topo_shape, Np,  Np_ovlp, &
                                   rank, comm_group, randField)
        implicit none
        !INPUT
        integer, dimension(3), intent(in) :: topo_pos, topo_shape, Np, Np_ovlp
        integer, intent(in) :: rank, comm_group
        !INPUT OUTPUT 
        double precision, dimension(:,:,:), intent(inout) :: randField
        !LOCAL
        double precision :: avg, std_dev
        integer, dimension(3) :: maxPos
        double precision :: sumRF, sumRFsquare
        double precision :: totalSumRF, totalSumRFsquare
        integer(kind = 8) :: xNTotal
        integer :: code

        maxPos = Np
        where(topo_pos /= (topo_shape-1)) maxPos = Np - Np_ovlp
        xNTotal = product((topo_shape*(Np-Np_ovlp)) +Np_ovlp)
        print*, "xNTotal = ", xNTotal

        !AVERAGE
        sumRF = sum(randField(1:maxPos(1),1:maxPos(2),1:maxPos(3))) 
        call MPI_ALLREDUCE (sumRF,totalSumRF,1,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,comm_group,code)
        avg = sumRF/dble(xNTotal)
        if(rank==0) print*, "avg BEFORE = ", avg
        randField = randField - avg

        !VARIANCE
        sumRFsquare = sum(randField(1:maxPos(1),1:maxPos(2),1:maxPos(3))**2d0)
        call MPI_ALLREDUCE (sumRFsquare,totalSumRFsquare,1,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,comm_group,code)
        std_dev = sqrt(sumRFsquare/dble(xNTotal))
        if(rank==0) print*, "std_dev BEFORE = ", std_dev
        randField = randField/std_dev
                
        end subroutine normalize_field

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
        nOvlpMax = 1 !We put 1 for security, it could probably be 0
        do dir = 1, size(neigh_rank)
            if(neigh_rank(dir) == rank) cycle 
            minP(:) = 1
            maxP(:) = Np
            where(neigh_shift(:,dir) == 1 ) minP = Np - Np_ovlp + 1
            where(neigh_shift(:,dir) == -1) maxP = Np_ovlp
            if(all(maxP > minP)) nOvlpMax = nOvlpMax + product(maxP - minP +1)
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
            if(any(maxP <= minP)) snd = .false. !Degenerated cases

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
            if(any(maxP <= minP)) rcv = .false. !Degenerated cases

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

                if(Np_ovlp(ii)<2) cycle
    
                if(topo_pos(ii) + shift < 0 &
                   .or. topo_pos(ii) + shift > (topo_shape(ii) - 1)) cycle 
    
                if(allocated(pattern)) deallocate(pattern)
                allocate(pattern(Np_ovlp(ii)))
                pattern = [(dble(kk), kk=0, Np_ovlp(ii)-1)]/dble((Np_ovlp(ii)-1))
                if(jj == 1) pattern = pattern(size(pattern):1:-1) !Pattern inversion
                !Chosing type of partition of unity
                if(partitionType == 1) then
                    pattern = ((1.0D0 + cos(PI*(pattern)))&
                                / 2.0D0)
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

    !--------------------------------------------------
    !--------------------------------------------------
    !--------------------------------------------------
    !--------------------------------------------------
    subroutine calculate_random_nb(nb)

        implicit none
        !OUTPUT
        integer, intent(out) :: nb
        !LOCAL
        integer :: i
        integer :: n
        integer :: clock, tempClock

        call system_clock(COUNT=clock)
        tempClock = clock
        do while (clock == tempClock)
            call system_clock(COUNT=clock)
        end do

        nb = clock + 37

    end subroutine calculate_random_nb
    !---------------------------------------
    !---------------------------------------
    !---------------------------------------
    !---------------------------------------
    subroutine init_random_seed(seedIn, seedStart)
        implicit none
        !INPUT
        integer, dimension(1:), optional, intent(in) :: seedIn
        integer, optional, intent(in) :: seedStart
        !LOCAL
        integer, dimension(:), allocatable :: seed

        if(present(seedIn)) then
            call random_seed(PUT = seedIn)
        else if (present(seedStart)) then
            call calculate_random_seed(seed, seedStart)
            call random_seed(PUT = seed)
            deallocate(seed)
        else
            call calculate_random_seed(seed)
            call random_seed(PUT = seed)
            deallocate(seed)
        end if

    end subroutine init_random_seed
    !------------------------------------------------
    !------------------------------------------------
    !------------------------------------------------
    !------------------------------------------------
    subroutine multiVariateTransformation (avg, std_dev, margiFirst, &
                                           randField)

        implicit none

        !INPUT
        integer          , intent(in) :: margiFirst
        double precision , intent(in) :: avg, std_dev

        !OUTPUT (IN)
        double precision, dimension(:,:,:), intent(inout) :: randField

        !LOCAL VARIABLES
        double precision :: normalVar, normalAvg
        integer          :: error, code

        select case (margiFirst)
        case(fom_GAUSSIAN)
            normalVar = (std_dev)**2d0
            normalAvg = avg
        case(fom_LOGNORMAL)
            if(avg <= 0.0D0) then
                write(*,*) ""
                write(*,*) "ERROR - when using lognormal fieldAvg should be a positive number"
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            normalVar = log(1d0 + (std_dev/avg)**2d0)
            normalAvg = log(avg) - normalVar/2d0
        case default
            print*, "margiFirst = ", margiFirst
            stop("First-order marginal not detected")
        end select

        randField = randField * sqrt(normalVar) &
                         + normalAvg;

        if (margiFirst == fom_LOGNORMAL) then
            randField = exp(randField)
        end if

    end subroutine multiVariateTransformation

end program main_ScaRL
