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
    !integer(kind = 8) :: pointsPerBlockIdeal = 250*250*250
    integer(kind = 8) :: pointsPerBlockIdeal = 1*1*1
       
    !LOCAL VARIABLES
    integer, dimension(3) :: Np, Np_temp
    integer :: comm_group
    integer, dimension(3) :: L
    integer, dimension(3) :: Np_ovlp, topo_shape, topo_pos
    double precision, dimension(3) :: xStep
    integer :: color, nb_procs_tmp, rank_tmp, code
    character (len=1024) :: command
    double precision, dimension(10) :: time_trace
    character (len=20), dimension(10) :: time_label
    integer :: time_count
    integer :: s
    double precision :: time_tic, time_toc, time_init

    integer :: gen_meth
    logical :: c_info_file, d_sample, one_file, one_hdf_ds
    character (len=1024) :: input_file
    integer, dimension(3) :: topo_shape_in
    logical ::extra_config
    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, comm_group, rank, nb_procs)

    !time_label(time_count) = "Initial"
    time_init = MPI_Wtime() 
    time_tic = time_init

    !print*, "TESTE FFTw";
    !call test_FFTw()
    call read_config_ScaRL("./config_ScaRL.txt", gen_meth, &
                          c_info_file, d_sample, one_file, one_hdf_ds, &
                          input_file, rank, comm_group) 

    call read_input_ScaRL(input_file, rank, &
                          nSamples, output_folder, &
                          xMinGlob, xMaxGlob, &
                          corrL, corrMod, margiFirst, &
                          seedBase, &
                          output_name, &
                          avg, std_dev, overlap, &
                          pointsPerCorrL, comm_group)
    
    call read_extra_config_ScaRL("extra_config_ScaRL.txt", extra_config, &
                                 topo_shape_in, rank, comm_group)
    
    !Creating Output Folder 
    if(rank == 0) then
        command = 'mkdir -pv '// trim(adjustL(output_folder)) 
        call system(command)
    end if
    call MPI_BARRIER (comm_group ,code)

    do s = 1, nSamples !----------------------------
    time_toc = MPI_Wtime()
    time_count = 1
    time_label(time_count) = "Read_input"
    time_trace(time_count) = time_toc-time_tic
    time_tic = MPI_Wtime()
   
    if(rank == 0) print*, "==================================================== "
    if(rank == 0) print*, "Generating SAMPLE --",s 
    if(rank == 0) print*, " "
    if(rank == 0) print*, " "
    if(rank == 0) print*, " "
     
    call verify_inputs(xMinGlob(:,s), xMaxGlob(:,s), corrL(:,s), & 
                       overlap(:,s), avg(s), std_dev(s), corrMod(s), margiFirst(s),&
                       pointsPerCorrL(:,s), pointsPerBlockIdeal, gen_meth)

    !Defining L, Np_ovlp and xStep (Np stands for Number of points)
    xStep   = corrL(:,s)/(dble(pointsPerCorrL(:,s)-1))  
    !print*, "xStep = ", xStep
    L       = 1+ceiling((xMaxGlob(:,s)-xMinGlob(:,s))/xStep)
    !print*, "xMaxGlob(:,s)-xMinGlob(:,s) = ", xMaxGlob(:,s)-xMinGlob(:,s)
    !print*, "ceiling = ", ceiling((xMaxGlob(:,s)-xMinGlob(:,s))/corrL(:,s))
    !print*, "L = ", L
    if(gen_meth == FFT_MPI) then
        if(rank == 0) print*, "WARNING!! In FFT MPI overlap is ignored"
        Np_ovlp = 0
    else
        Np_ovlp = 1+ceiling(overlap(:,s)*dble(pointsPerCorrL(:,s)-1))
    end if
    !print *, "rank, Np_ovlp =", rank, Np_ovlp
    !print *, "rank, overlap =", rank, overlap(:,s)
    !print *, "rank, pointspcorrl =", rank, pointsPerCorrL(:,s)
    !print *, "xMinGlob   =", L
    !print *, "xMaxGlob   =", L
    !print *, "L          =", L
    !Finding place in topology and defining the number of points per processor (Np)
    if(extra_config) then
        if( product(topo_shape_in) /= nb_procs ) stop 'ERROR! extra_config_ScaRL.txt found but procs grid dimensions do not agree with the total number of procs in MPI process'
        topo_shape = topo_shape_in 
    else 
        call decide_topo_shape(nb_procs, L, Np_ovlp, &
                           pointsPerBlockIdeal, & 
                           topo_shape, rank, corrMod(s),MPI_COMM_WORLD, &
                           gen_meth)
    end if
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
        !print *, "L   rg", rank," =", L
        !print *, "topo_shape   rg", rank," =", topo_shape
        !print *, "Np_ovlp  rg", rank," =", Np_ovlp

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

        !print *, "Np  rg", rank," =", Np
        L = Np*topo_shape - ((topo_shape-1)*Np_ovlp)
        if(rank == 0) print *, "nb_procs   =", nb_procs
        if(rank == 0) print *, "topo_shape =", topo_shape
        if(rank == 0) print *, "L          =", L
        if(rank == 0) print *, "Np         =", Np
        if(rank == 0) print *, "Np_ovlp    =", Np_ovlp
        if(rank == 0) print *, "product(Np)=", product(Np)
        if(rank == 0) print *, "xStep      =", xStep

        call get_topo_pos(rank, topo_shape, topo_pos)
        
        time_toc = MPI_Wtime()
        time_count = time_count + 1
        time_label(time_count) = "Topology"
        time_trace(time_count) = time_toc-time_tic
        time_tic = MPI_Wtime()

       !Creating Fields
        call create_fields(Np, Np_ovlp, L, pointsPerCorrL(:,s), rank, &
                           nb_procs, topo_pos, topo_shape, &
                           xStep, xMinGlob(:,s), &
                           corrL(:,s), corrMod(s), &
                           seedBase(s), & 
                           margiFirst(s), avg(s), std_dev(s), &
                           comm_group, &
                           output_name(s), &
                           output_folder, &
                           time_count, time_label, time_trace, time_init, &
                           gen_meth, &
                           c_info_file, d_sample, one_file, one_hdf_ds)
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
                         pointsPerCorrL, pointsPerBlockIdeal, gen_meth)

implicit none 
double precision, dimension (3), intent(in) :: xMinGlob
double precision, dimension (3), intent(in) :: xMaxGlob
double precision, dimension (3), intent(in) :: corrL
double precision, dimension (3), intent(in) :: overlap
double precision, intent(in) :: avg
double precision, intent(in) :: std_dev
integer, intent(in) :: corrMod, margiFirst, gen_meth
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

if((gen_meth /= FFT) .and. &
    (gen_meth /= FFT_MPI)) then
    print*, "gen_meth = ", gen_meth
    stop("ERROR: Generation Method  not implemented")
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

!if(pointsPerBlockIdeal < 8) then
!    print*, "pointsPerBlockIdeal = ", pointsPerBlockIdeal
!    stop("ERROR: pointsPerBlockIdeal < 8")
!end if 

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
                         output_name_in, res_folder, &
                         time_count, time_label, time_trace, time_init, &
                         gen_meth, &
                         c_info_file, d_sample, one_file, one_hdf_ds)
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
    character(len=*), intent(in) :: output_name_in, res_folder
    double precision, intent(in) :: time_init
    integer, intent(in) :: gen_meth
    logical, intent(in) :: c_info_file, d_sample, one_file, one_hdf_ds

    !OUTPUT
    integer, intent(inout) :: time_count 
    character(len=*), dimension(:), intent(inout) :: time_label
    double precision, dimension(:), intent(inout) :: time_trace
     
    !LOCAL
    double precision, dimension(:, :, :), allocatable :: k_mtx
    character(len=1024) :: output_name, HDF5_name, XMF_name
    double precision, dimension(3) :: coord_0, coord_N
    integer         , dimension(3) :: pos_0, pos_N
    double precision, dimension(3) :: xMaxGlob, xRangeGlob
    double precision, dimension(3,nb_procs) :: coord_0_list
    character(len=1024), dimension(nb_procs) :: HDF5_list 
    integer :: temp_rank
    integer, dimension(3) :: temp_origin
    integer, dimension(3) :: temp_topo_pos
    integer :: partition_type = 1
    integer :: seedStart
    double precision, dimension(3) :: xRange, overlap
    integer :: code
    double precision :: time_tic, time_toc
    character(len=10), dimension(3) :: strings
    integer, dimension(8) :: date_time
    character(len=50) :: date_str
    character(len=1024) :: filePath, info_file_name
    logical :: fileExists
    double precision :: time_final
    integer :: ii


    allocate(k_mtx(Np(1), Np(2), Np(3)))
    !print*, "L= ", L
    time_tic = MPI_Wtime()
    
    output_name = output_name_in
    
    if(one_file) then 
        HDF5_name = str_cat(output_name,".h5")
        XMF_name  = str_cat(output_name,".xmf")
    else
        if(rank == 0) print*, "HDF5 - One File Per Proc and One Dataset"
        HDF5_name = str_cat(output_name,"_proc_",trim(num2str(rank,5)),".h5")
        XMF_name  = str_cat(output_name,"_proc_",trim(num2str(rank,5)),".xmf")
    end if
    
    filePath = str_cat(res_folder,"/",HDF5_name)
    if(rank == 0) print*, "filePath = ", trim(filePath)
    inquire(file=filePath, exist=fileExists)
    if(rank == 0) print*, "     exists? ",fileExists
    
    if(fileExists) then
        if(rank == 0) then
            call date_and_time(strings(1), strings(2), strings(3), date_time)
            date_str = strings(1)(3:8)//"_"//strings(2)(1:6)
        end if
        
        call MPI_BCAST (date_str, len(date_str), &
                        MPI_CHARACTER, 0, comm_group, code)

        output_name = str_cat(output_name_in,"-",date_str)
        
        if(one_file) then 
            HDF5_name = str_cat(output_name,".h5")
            XMF_name  = str_cat(output_name,".xmf")
        else
            if(rank == 0) print*, "HDF5 - One File Per Proc and One Dataset"
            HDF5_name = str_cat(output_name,"_proc_",trim(num2str(rank,5)),".h5")
            XMF_name  = str_cat(output_name,"_proc_",trim(num2str(rank,5)),".xmf")
        end if
        
        if(rank == 0) print*, "NEW output_name    = ", trim(output_name)
    
    else
        if(rank == 0) print*, "output_name    = ", trim(output_name)
    
    end if
     
    if(rank == 0) print*, "HDF5 output proc 0 = ", trim(HDF5_name)
    if(rank == 0) print*, "XMF  output proc 0 = ", trim(XMF_name)

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
    xRangeGlob = xMaxGlob - xMinGlob 

    !k_mtx(:,:,:) = dble(rank)
    k_mtx(:,:,:) = 1d0
!    if(rank == 0) print *, "Np                  = ", Np
!    if(rank == 0) print *, "Np_ovlp             = ", Np_ovlp
!    if(rank == 0) print *, "topo_pos            = ", topo_pos
!    if(rank == 0) print *, "pos_0               = ", pos_0
!    if(rank == 0) print *, "pos_N               = ", pos_N
!    if(rank == 0) print *, "coord_0             = ", coord_0
!    if(rank == 0) print *, "coord_N             = ", coord_N
!    if(rank == 0) print *, "xMinGlob            = ", xMinGlob
!    if(rank == 0) print *, "xMaxGlob            = ", xMaxGlob
!    if(rank == 0) print *, "xRange (processor)  = ", xRange
!    if(rank == 0) print *, "xRange (total)      = ", xMaxGlob - xMinGlob
    do ii=0, nb_procs-1
        if(rank == ii) print *, "RANK ", rank, " ----------------"
        if(rank == ii) print *,  "Np                  = ", Np
        if(rank == ii) print *,  "Np_ovlp             = ", Np_ovlp
        if(rank == ii) print *,  "topo_pos            = ", topo_pos
        if(rank == ii) print *,  "pos_0               = ", pos_0
        if(rank == ii) print *,  "pos_N               = ", pos_N
        if(rank == ii) print *,  "coord_0             = ", coord_0
        if(rank == ii) print *,  "coord_N             = ", coord_N
        if(rank == ii) print *,  "xMinGlob            = ", xMinGlob
        if(rank == ii) print *,  "xMaxGlob            = ", xMaxGlob
        if(rank == ii) print *,  "xRange (processor)  = ", xRange
        if(rank == ii) print *,  "xRange (total)      = ", xMaxGlob - xMinGlob
        call MPI_BARRIER(comm_group, code)
    end do
    !print*, "2L= ", L
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
    
    if (gen_meth == FFT) then
        if(rank == 0) print*, "Generating sample using FFT  with localization" 
        call gen_std_gauss_FFT(k_mtx, Np, &
                               xRange, corrL, &
                               corrMod, &
                               seedStart, rank)
    else if (gen_meth == FFT_MPI) then
        if(rank == 0) print*, "Generating sample using FFT (only MPI) " 
        call gen_Std_Gauss_FFT_MPI(k_mtx, &
                                   xRangeGlob, corrL, corrMod, &
                                   seedStart, rank, L, Np, comm_group)
    end if
    !call gen_std_gauss_Shino_FFT(k_mtx, Np, &
    !                       xRange, corrL, &
    !                       pointsPerCorrL, corrMod, &
    !                       seedStart, rank)
    
    time_toc = MPI_Wtime()
    time_count = time_count + 1
    time_label(time_count) = "Generation"
    time_trace(time_count) = time_toc-time_tic
    time_tic = MPI_Wtime()

    if (gen_meth == FFT) then
        if(rank == 0) print*, "apply_UnityPartition " 
        call apply_UnityPartition_mtx(Np, Np_ovlp,&
                                      partition_type, &
                                      topo_pos, topo_shape, &
                                      k_mtx)
    else if (gen_meth == FFT_MPI) then
        if(rank == 0) print*, "FFT_MPI. No need for Unity Partition. " 
    end if
    
    time_toc = MPI_Wtime()
    time_count = time_count + 1
    time_label(time_count) = "Apply_unit"
    time_trace(time_count) = time_toc-time_tic
    time_tic = MPI_Wtime()
   
    if (gen_meth == FFT) then
        if(rank == 0) print*, "add_overlap " 
        call add_overlap(k_mtx, Np, Np_ovlp, rank, &
                       nb_procs, topo_pos, topo_shape, &
                       comm_group)
    else if (gen_meth == FFT_MPI) then
        if(rank == 0) print*, "FFT_MPI. No overlap to add. " 
    end if
    
    time_toc = MPI_Wtime()
    time_count = time_count + 1
    time_label(time_count) = "Add_overlap"
    time_trace(time_count) = time_toc-time_tic
    call MPI_Barrier(comm_group,code)
    time_tic = MPI_Wtime()
  
    !if(rank == 0) print*, "normalize_field " 
    !call normalize_field(topo_pos, topo_shape, Np, Np_ovlp, &
    !                    rank, comm_group, k_mtx)

    time_toc = MPI_Wtime()
    time_count = time_count + 1
    time_label(time_count) = "Normalize"
    time_trace(time_count) = time_toc-time_tic
    time_tic = MPI_Wtime()
   
   if(rank == 0) print*, "multivariateTransformation " 
   call multiVariateTransformation(avg, std_dev, margiFirst, &
                                   k_mtx)

    time_toc = MPI_Wtime()
    time_count = time_count + 1
    time_label(time_count) = "Multivar_trans"
    time_trace(time_count) = time_toc-time_tic
    time_tic = MPI_Wtime()
   
   if(rank == 0) print*, "maxval(k_mtx) AFTER = ", maxval(k_mtx) 
   if(rank == 0) print*, "minval(k_mtx) AFTER = ", minval(k_mtx)

   if(rank == 0) print*, "write HDF5 "


   if(one_file) then 
       !HDF5_name = str_cat(output_name,".h5")
       !XMF_name  = str_cat(output_name,".xmf")
       if(one_hdf_ds) then
           if(rank == 0) print*, "HDF5 - One File and One Dataset"
           !print*, "1 L= ", L
           call write_hdf5_multi_proc_3D_1ds(coord_0, &
                            coord_N,     &
                            xMinGlob, xMaxGlob, &
                            pos_0, pos_N, &
                            L, Np, Np_ovlp, &
                            k_mtx,       &
                            HDF5_name,  &
                            XMF_name,  &
                            res_folder,           &
                            rank, nb_procs,           &
                            comm_group)
       else 
           if(rank == 0) print*, "HDF5 - One File and Several Datasets"
           call write_hdf5_multi_proc_3D(coord_0, &
                            coord_N,     &
                            xMinGlob, xMaxGlob, &
                            pos_0, pos_N, &
                            k_mtx,       &
                            HDF5_name,  &
                            res_folder,           &
                            rank, nb_procs,           &
                            comm_group)
       end if
       if(rank == 0) then
           call write_HDF5_attributes(str_cat(res_folder,"/",output_name,".h5"), &
                        nb_procs, 3, 1, FFT, seedBase, &
                        corrMod, margiFirst, &
                        topo_shape, &
                        xMinGlob, xMaxGlob, xStep, corrL, overlap, &
                        L, Np, Np_ovlp, .false.)
       end if
   else
       if(rank == 0) print*, "HDF5 - One File Per Proc and One Dataset"
       !HDF5_name = str_cat(output_name,"_proc_",trim(num2str(rank,5)),".h5")
       !XMF_name  = str_cat(output_name,"_proc_",trim(num2str(rank,5)),".xmf")
       call write_hdf5_single_proc_3D(coord_0, &
                                  coord_N, &
                                  xMinGlob, xMaxGlob, &
                                  pos_0, pos_N, &
                                  k_mtx, &
                                  HDF5_name, &
                                  XMF_name, &
                                  res_folder, &
                                  rank)
       call write_HDF5_attributes(str_cat(res_folder,"/",HDF5_name), &
                    nb_procs, 3, 1, FFT, seedBase, &
                    corrMod, margiFirst, &
                    topo_shape, &
                    xMinGlob, xMaxGlob, xStep, corrL, overlap, &
                    L, Np, Np_ovlp, .false.)
   end if

   if(rank==0 .and. (.not.(one_file .and. one_hdf_ds))) then
       if(rank == 0) print*, "WRITING GLOBAL XMF"
       print*, "HERE"
       do temp_rank = 0, nb_procs-1
           call get_topo_pos(temp_rank,  &
                             topo_shape, &
                             temp_topo_pos)
           HDF5_list(temp_rank+1) = &
            str_cat(output_name,"_proc_",trim(num2str(temp_rank,5)),".h5")
           if(one_file)  HDF5_list(temp_rank+1) = &
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
   
    time_toc = MPI_Wtime()
    time_count = time_count + 1
    time_label(time_count) = "Write_file"
    time_trace(time_count) = time_toc-time_tic
    time_tic = MPI_Wtime()

    call MPI_BARRIER(comm_group, code)
    time_final = MPI_Wtime()
    time_count = time_count + 1
    time_label(time_count) = "Total_Wall_time"
    time_trace(time_count) = time_final-time_init

   if(rank==0) print*, "Writing times on ", trim(HDF5_name)
   call write_time_on_hdf5(time_label(1:time_count), &
                           time_trace(1:time_count), &
                           HDF5_name, &
                           res_folder, &
                           rank, nb_procs, comm_group, one_file)

   if(c_info_file) then
       if(rank==0) print*, "Writing INFO file"

       if(rank==0) then
           info_file_name = str_cat("INFO-",output_name,".h5")
        
           filePath = str_cat(res_folder,"/",info_file_name)
           if(rank == 0) print*, "filePath = ", trim(filePath)
           inquire(file=filePath, exist=fileExists)
           if(rank == 0) print*, "     exists? ",fileExists
           
           if(fileExists) then
               call date_and_time(strings(1), strings(2), strings(3), date_time)
               date_str = strings(1)(3:8)//"_"//strings(2)(1:6)
               info_file_name = str_cat("INFO-",output_name,"-",date_str,".h5")
           end if

           print*, "info_file_name = ", trim(info_file_name)
           
           call write_info_file(info_file_name, res_folder, rank)
           print*, "Writing attributes"
           call write_HDF5_attributes(&
                        str_cat(res_folder,"/",info_file_name), &
                        nb_procs, 3, 1, FFT, seedBase, &
                        corrMod, margiFirst, &
                        topo_shape, &
                        xMinGlob, xMaxGlob, xStep, corrL, overlap, &
                        L, Np, Np_ovlp, .false.)
       end if
       
       call MPI_BCAST (info_file_name, len(info_file_name,), &
                       MPI_CHARACTER, 0, comm_group, code)

       if(rank==0) print*, "Writing times"
       call write_time_on_hdf5(time_label(1:time_count), &
                               time_trace(1:time_count), &
                               info_file_name, &
                               res_folder, &
                               rank, nb_procs, comm_group, &
                               oneFile=.true.)
   end if

   if(d_sample .and. rank == 0) then
      print*, "WARNING!! Deleting sample files (only INFO files will stand)"
      !call system("rm SAMPLES/S*")
      call system(str_cat("rm "//res_folder,"/",output_name,"*"))

   end if
   
   if(allocated(k_mtx)) deallocate(k_mtx)

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
                             topo_shape, rank, corrMod, comm, &
                             gen_meth)

    implicit none
    !INPUT
    integer, intent(in) :: nb_procs
    integer, dimension(3), intent(in) :: L, nPointsOvlp
    integer(kind = 8), intent(in) :: pointsPerBlockIdeal
    integer, intent(in) :: rank, corrMod, comm, gen_meth
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
    !FFTW_MPI
    integer(C_INTPTR_T) :: L1, L2, L3
    integer(C_INTPTR_T) :: local_LastDim
    integer(C_INTPTR_T) :: local_LD_offset
    integer(C_INTPTR_T) :: alloc_local
    integer :: count_proc, max_procs_dim3

    if(gen_meth == FFT_MPI) then
        L1 = L(1)
        L2 = L(2)
        L3 = L(3)
        alloc_local = fftw_mpi_local_size_3d(L3, L2, L1, comm, &
                                local_LastDim, local_LD_offset) !FOR MPI
        count_proc=1
        if(local_LastDim == 0) then
            count_proc=0
            ! "PROCESSOR IGNORED BY FFTW DECOMPOSITION")
        end if
                
        call MPI_ALLREDUCE (count_proc,max_procs_dim3,1,MPI_INTEGER, &
                            MPI_SUM,comm,code)
        topo_shape(:) = 1
        topo_shape(3) = max_procs_dim3
    else
        nPointsBase = L
        nBlocks = 1
        nFieldsOK = .false.
        do while (.not. nFieldsOK)
        
            nPoints = nPointsBase + (nBlocks - 1)*nPointsOvlp
            nPointsPerProc = ceiling(dble(nPoints)/dble(nBlocks))
            nTotalPointsProc = product(int(nPointsPerProc,8))
            if(nTotalPointsProc <= pointsPerBlockIdeal) then
                nFieldsOK = .true.
            else if(product(nBlocks) >= nb_procs) then 
                nFieldsOK = .true. 
            else
                pos = maxloc(nPointsPerProc, 1)
                nBlocks(pos) = nBlocks(pos) + 1
            end if

            !print*, "nPoints = ", nPoints
            !print*, "nPointsPerProc = ", nPointsPerProc
            !print*, "nBlocks = ", nBlocks
        
        end do

        nBlocksIdeal = nBlocks
        if(rank == 0) print*, "L = ", L
        if(rank == 0) print*, "nBlocksIdeal = ", nBlocksIdeal

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
    end if

    if(rank == 0) print*, "topo_shape = ", topo_shape

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
integer(kind = 8), dimension(3) :: topo_shape_8, Np_8, Np_ovlp_8
integer :: code

maxPos = Np
where(topo_pos /= (topo_shape-1)) maxPos = Np - Np_ovlp
topo_shape_8 = topo_shape
Np_8 = Np
Np_ovlp_8 = Np_ovlp
xNTotal = product((topo_shape_8*(Np_8-Np_ovlp_8))+Np_ovlp_8)
if(rank==0) print*, "    Normalizing average and stdDev of Gaussian ensemble"
if(rank==0) print*, "    xNTotal = ", xNTotal

!AVERAGE
sumRF = sum(randField(1:maxPos(1),1:maxPos(2),1:maxPos(3))) 
call MPI_ALLREDUCE (sumRF,totalSumRF,1,MPI_DOUBLE_PRECISION, &
                    MPI_SUM,comm_group,code)
avg = totalSumRF/dble(xNTotal)
randField = randField - avg

!VARIANCE
sumRFsquare = sum(randField(1:maxPos(1),1:maxPos(2),1:maxPos(3))**2d0)
call MPI_ALLREDUCE (sumRFsquare,totalSumRFsquare,1,MPI_DOUBLE_PRECISION, &
                    MPI_SUM,comm_group,code)
std_dev = sqrt(totalSumRFsquare/dble(xNTotal))
randField = randField/std_dev
print*, "    loc_sumRF       rank ", rank, "  = ", sumRF
print*, "    loc_sumRFsquare rank ", rank, "  = ", sumRFsquare
!if(rank==0) print*, "    loc_sumRF        = ", sumRF
!if(rank==0) print*, "    loc_sumRFsquare  = ", sumRFsquare
if(rank==0) print*, "    glob_sumRF       = ", totalsumRF
if(rank==0) print*, "    glob_sumRFsquare = ", totalsumRFsquare
if(rank==0) print*, "    avg BEFORE       = ", avg
if(rank==0) print*, "    std_dev BEFORE   = ", std_dev
if(rank==0) print*, "    avg AFTER (supposed) = ", 0d0
if(rank==0) print*, "    std_dev AFTER (supposed) = ", 1d0
        
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
    double precision, dimension(:, :, :), allocatable :: RF_temp1, RF_temp2
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
    integer(kind=8):: double_size
    integer :: double_size_short
    integer(kind=8) :: totalSize, overHead, overEst, bufferSize
    double precision, dimension(:), allocatable :: buffer
    integer :: request
    integer, dimension(MPI_STATUS_SIZE) :: statut
    integer :: tag
    logical :: snd, rcv
    integer :: BufDT_size
    integer :: ii, jj, kk, cnt
    integer :: rang_test=-1,r_sender = 1, r_receiver = 6

    allocate(RF_temp1(Np(1), Np(2), Np(3)))
    allocate(RF_temp2(Np(1), Np(2), Np(3)))

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
                
                if(neigh_rank(cnt) /= rank .or. op_neigh_rank(cnt) /= rank)then 
                if(rang_test == rank) print*, "rk ", rank, "neigh shift", neigh_shift(:,cnt), "---------------"
                if(rang_test == rank .and. neigh_rank(cnt) /= rank) &
                      print*, "neigh = ", neigh_rank(cnt)
                if(rang_test == rank .and. op_neigh_rank(cnt) /= rank) &
                      print*, "op_ng = ", op_neigh_rank(cnt) 
                if(rang_test == rank) print*, "  " 
                end if

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
        if(all(maxP > minP)) nOvlpMax = nOvlpMax + product(int(maxP,8) - int(minP,8) +1)
    end do

    !Buffer allocation
    overEst = 2
    call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,double_size_short,code)
    double_size = double_size_short
    overHead = 1+int(MPI_BSEND_OVERHEAD,8)/double_size
    bufferSize = overEst*(nOvlpMax+overHead)
    if(int(double_size*bufferSize) < 0) then
        print*, "ERROR!!!!  Overflow of buffer_MPI (>2147483647, for INT_32)"
        print*, "    Points per processor      = ", Np
        print*, "    Overlap points            = ", Np_ovlp
        print*, "    Total points in overlap   = ", nOvlpMax
        print*, "    overHead                  = ", overHead
        print*, "    buffer_MPI                = ", bufferSize*double_size
        print*, "    buffer_MPI/max_buffer_MPI = ", 100d0*dble(bufferSize*double_size)/(2e9),"%"
        stop ">>> Consider reducing the overlaping zone or using more processors"
    end if
    allocate(buffer(bufferSize))
    call MPI_BUFFER_ATTACH(buffer, int(double_size*bufferSize),code)
    if(code /= 0) then
        print*, "ERROR in MPI_BUFFER_ATTACH"
        print*, "MPI_BUFFER_ATTACH code = ",code
        stop(" ")
    end if
    !Communications
    RF_temp1 = 0d0

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
        !if(rank == r_sender) print*, "Np = ", Np
        !print*, "Np = ", Np
        !if(rank == r_sender) print*, "Np_ovlp = ", Np_ovlp
        !if(rank == r_sender) print*, "direction = ", dir
        !if(rank == r_sender) print*, "dirShift  = ", dirShift
        !if(rank == r_sender) print*, "minP  = ", minP
        !if(rank == r_sender) print*, "maxP  = ", maxP
        
        !if(rank == r_sender) then
        !    print*, "ANALYSE rank = ", neigh_rank(dir)
        !end if

        !if(rank == r_sender) print*, "snd 1= ", snd
        if(neigh_rank(dir) == rank) snd = .false. !Check if this direction exists
        !if(rank == r_sender) print*, "snd 2= ", snd
        if(any(maxP <= minP)) snd = .false. !Degenerated cases
        !if(rank == r_sender) print*, "snd 3= ", snd

        !TEST
        !if(rank /= r_sender) snd=.false.
        !if(rank == r_sender) print*, "snd 4= ", snd
        !if(neigh_rank(dir) /= r_receiver) snd=.false.
        !if(rank == r_sender) print*, "snd 5= ", snd

        if(snd) then
            !if(rank == r_sender) print*, " SENDING =============="
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

        !TEST
        !if(rank /= r_receiver) rcv=.false.
        !if(neigh_rank(op_dir) /= r_sender) rcv=.false.

        if(rcv) then
            !if(rank == r_receiver) print*, " RECEIVING =============="

            tag  = rank
            RF_temp2 = 0.0D0

            call MPI_RECV (RF_temp2(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3)), &
                int(totalSize), MPI_DOUBLE_PRECISION, &
                neighRank, tag, comm_group, statut, code)
                
            RF_temp1(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3)) = &
                RF_temp1(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))   &
                + RF_temp2(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))
            
            if(rank == rang_test) print*, "  RECVER: rank ", rank
            if(rank == rang_test) print*, "  FROM rank ", neighRank
            if(rank == rang_test) print*, "  minP = ", minP
            if(rank == rang_test) print*, "  maxP = ", maxP
            if(rank == rang_test) print*, "  op_dir_shft = ", -dirShift
            if(rank == rang_test) print*, "  totalSize = ", totalSize
            if(rank == rang_test) print*, "  tag = ", tag
            if(rank == rang_test) print*, "  CONTENT = ", RF_temp2(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))
        else
        
        end if
    end do

    RF = RF + RF_temp1

    BufDT_size = int(double_size*bufferSize)
    call MPI_BUFFER_DETACH (buffer,BufDT_size,code)
    
    if(allocated(buffer)) deallocate(buffer)
    if(allocated(RF_temp1)) deallocate(RF_temp1)
    if(allocated(RF_temp2)) deallocate(RF_temp2)

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
    double precision, dimension(:, :, :), allocatable :: unityPartition
    integer :: shift

    allocate(unityPartition(Np(1), Np(2), Np(3)))
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
    if(allocated(unityPartition)) deallocate(unityPartition)

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
