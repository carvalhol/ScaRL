module read_input

    use str_functions
    use mpi
    implicit none
    

    type :: property_RF
        integer :: mat
        character(len=1024) :: name
        double precision :: avg
        integer :: corrMod
        double precision, dimension(3) :: corrL
        integer :: margiF
        double precision :: CV
        integer :: seedStart
        double precision, dimension(3) :: bbox_min, bbox_max

    end type property_RF

contains
        !--------------------------------------------------------------
        !--------------------------------------------------------------
        !--------------------------------------------------------------
        !--------------------------------------------------------------
        subroutine read_input_ScaRL(file_path_in, rank, &
                                    nSamples, output_folder, &
                                    xMinGlob, xMaxGlob, &
                                    corrL, corrMod, margiFirst, &
                                    seedBase, &
                                    output_name, &
                                    avg, std_dev, overlap, &
                                    pointsPerCorrL, comm_group)
            implicit none
            !INPUT
            character(len=*), intent(in) :: file_path_in
            integer, intent(in) :: rank
            !OUTPUT
            integer, intent(out) :: nSamples
            character (len=1024) , intent(out):: output_folder
            double precision, dimension (:,:), allocatable, intent(out) :: xMinGlob
            double precision, dimension (:,:), allocatable, intent(out) :: xMaxGlob
            double precision, dimension (:,:), allocatable, intent(out) :: corrL
            integer, dimension(:), allocatable, intent(out) :: corrMod, margiFirst
            integer, dimension(:), allocatable, intent(out) :: seedBase
            character (len=1024), dimension(:), allocatable, intent(out) :: output_name
            double precision, dimension(:), allocatable, intent(out) :: avg, std_dev
            double precision, dimension (:,:), allocatable, intent(out) :: overlap
            integer, dimension(:,:), allocatable, intent(out) :: pointsPerCorrL
            integer, intent(in) :: comm_group
            !LOCAL
            character(len=1024) :: buffer,cur_folder,file_path
            integer :: fid
            integer :: s, error, code
        integer :: char_size
        integer(kind=8) :: bufferSize
        character(len=2048) :: bufferChar 
        

        if(rank == 0) then
            call getcwd(cur_folder)
            print*, "Working on directory: ", trim(adjustL(cur_folder))
            file_path = str_cat(cur_folder,"/random.spec")

            if(fileExist(file_path)) then

                if(rank == 0) print*, "WARNING: 'random.spec' found, generating files for SEM"   
                file_path = str_cat(cur_folder,"/domains.txt")
                if(.not. fileExist(file_path)) then
                    print*, "'domains.txt'NOT FOUND (Did you run the mesher?"
                    call MPI_ABORT(comm_group, error, code)
                end if
                call read_input_SEM(cur_folder, rank, &
                                    nSamples, output_folder, &
                                    xMinGlob, xMaxGlob, &
                                    corrL, corrMod, margiFirst, &
                                    seedBase, &
                                    output_name, &
                                    avg, std_dev, overlap, &
                                    pointsPerCorrL)
            else
                file_path = str_cat(cur_folder,"/",file_path_in)

                if(rank == 0) print*, "READING INPUT FROM ", trim(adjustL(file_path))   

                if(.not. fileExist(file_path)) then
                    print*, "'input_Scarl.txt' NOT FOUND"
                    call MPI_ABORT(comm_group, error, code)
                end if
                
                    print*, "BEFORE"
                fid = 50 
                open (unit = fid , file = file_path, action = 'read', iostat=error)
                !open (unit = fid , file = "./input_ScaRL.txt", action = 'read', iostat=error)
                   
                    print*, "error = ", error 
                   
                    buffer = getLine(fid, '#') !output folder
                    read(buffer,*) output_folder
                    if(rank == 0) print*, "output folder = ", trim(adjustL(output_folder))

                    buffer = getLine(fid, '#') !Number of Samples
                    read(buffer,*) nSamples
                    if(rank == 0) print*, "Number of Samples = ", nSamples

                    allocate(output_name(nSamples))
                    allocate(xMinGlob(3,nSamples))    
                    allocate(xMaxGlob(3,nSamples))
                    allocate(corrL(3,nSamples))
                    allocate(overlap(3,nSamples))
                    allocate(pointsPerCorrL(3,nSamples))
                    allocate(corrMod(nSamples))
                    allocate(margiFirst(nSamples))
                    allocate(avg(nSamples))
                    allocate(std_dev(nSamples))
                    allocate(seedBase(nSamples))

                    do s = 1, nSamples
                    
                    if(rank == 0) print*, " "
                    if(rank == 0) print*, " "
                    if(rank == 0) print*, "SAMPLE ", s, "----------"
                    
                    buffer = getLine(fid, '#') !sample name
                    read(buffer,*) output_name(s)
                    if(rank == 0) print*, "sample name = ", trim(adjustL(output_name(s)))
                    
                    buffer = getLine(fid, '#') !xMinGlob
                    read(buffer,*) xMinGlob(:,s)
                    if(rank == 0) print*, "xMinGlob = ", xMinGlob(:,s)

                    buffer = getLine(fid, '#') !xMaxGlob
                    read(buffer,*) xMaxGlob(:,s)
                    if(rank == 0) print*, "xMaxGlob = ", xMaxGlob(:,s)
                    
                    buffer = getLine(fid, '#') !corrL
                    read(buffer,*) corrL(:,s)
                    if(rank == 0) print*, "corrL = ", corrL(:,s)
                    
                    buffer = getLine(fid, '#') !overlap
                    read(buffer,*) overlap(:,s)
                    if(rank == 0) print*, "overlap = ", overlap(:,s)
                    
                    buffer = getLine(fid, '#') !pointsPerCorrL
                    read(buffer,*) pointsPerCorrL(:,s)
                    if(rank == 0) print*, "pointsPerCorrL = ", pointsPerCorrL(:,s)
                    
                    buffer = getLine(fid, '#') !corrMod
                    read(buffer,*) corrMod(s)
                    if(rank == 0) print*, "corrMod = ", corrMod(s)
                    
                    buffer = getLine(fid, '#') !margiFirst
                    read(buffer,*) margiFirst(s)
                    if(rank == 0) print*, "first order marginal = ", margiFirst(s)
                    
                    buffer = getLine(fid, '#') !avg
                    read(buffer,*) avg(s)
                    if(rank == 0) print*, "avg = ", avg(s)
                    
                    buffer = getLine(fid, '#') !standard deviation
                    read(buffer,*) std_dev(s)
                    if(rank == 0) print*, "std_dev = ", std_dev(s)
                    
                    buffer = getLine(fid, '#') !random seed
                    read(buffer,*) seedBase(s)
                    if(rank == 0) print*, "random seed = ", seedBase(s)
                    
                    end do

                close (fid)

            end if
        end if
        
        !Passing information to others
        call MPI_BCAST (nSamples, 1, MPI_INTEGER, 0, comm_group, code)
        !print*, "rank ", rank, "nSamples", nSamples   
        if(rank /= 0) then     
            allocate(output_name(nSamples))
            allocate(xMinGlob(3,nSamples))    
            allocate(xMaxGlob(3,nSamples))
            allocate(corrL(3,nSamples))
            allocate(overlap(3,nSamples))
            allocate(pointsPerCorrL(3,nSamples))
            allocate(corrMod(nSamples))
            allocate(margiFirst(nSamples))
            allocate(avg(nSamples))
            allocate(std_dev(nSamples))
            allocate(seedBase(nSamples))
        end if

      !  call MPI_TYPE_SIZE(MPI_CHARACTER,char_size,code)
      !  bufferSize = len(bufferChar)
      !  call MPI_BUFFER_ATTACH(bufferChar, int(char_size*bufferSize),code)
      !  if(code /= 0) then
      !      print*, "ERROR in MPI_BUFFER_ATTACH"
      !      print*, "MPI_BUFFER_ATTACH code = ",code
      !      stop(" ")
      !  end if

      !  print*, "HERE BCAST"
        call MPI_BCAST (output_name, len(output_name), MPI_CHARACTER, 0, comm_group, code)
        print*, "rank ", rank, "output_name(1): ", trim(output_name(1)) 
        call MPI_BCAST (xMinGlob, 3*nSamples, MPI_DOUBLE_PRECISION, 0, comm_group, code)
        call MPI_BCAST (xMaxGlob, 3*nSamples, MPI_DOUBLE_PRECISION, 0, comm_group, code)
        call MPI_BCAST (corrL,    3*nSamples, MPI_DOUBLE_PRECISION, 0, comm_group, code)
        call MPI_BCAST (overlap,  3*nSamples, MPI_DOUBLE_PRECISION, 0, comm_group, code)
        call MPI_BCAST (pointsPerCorrL, 3*nSamples, MPI_INTEGER   , 0, comm_group, code)
        call MPI_BCAST (corrMod   , nSamples, MPI_INTEGER   , 0, comm_group, code)
        call MPI_BCAST (margiFirst, nSamples, MPI_INTEGER   , 0, comm_group, code)
        call MPI_BCAST (avg   , nSamples, MPI_DOUBLE_PRECISION, 0, comm_group, code)
        call MPI_BCAST (std_dev   , nSamples, MPI_DOUBLE_PRECISION, 0, comm_group, code)
        call MPI_BCAST (seedBase  , nSamples, MPI_INTEGER   , 0, comm_group, code)
        
        end subroutine read_input_ScaRL
        
        !--------------------------------------------------------------
        !--------------------------------------------------------------
        !--------------------------------------------------------------
        !--------------------------------------------------------------
        subroutine read_input_SEM( cur_folder, rank, &
                                   nSamples, output_folder, &
                                   xMinGlob, xMaxGlob, &
                                   corrL, corrMod, margiFirst, &
                                   seedBase, &
                                   output_name, &
                                   avg, std_dev, overlap, &
                                   pointsPerCorrL)
            implicit none
            !INPUT
            character(len=*) :: cur_folder
            integer, intent(in) :: rank
            !OUTPUT
            integer, intent(out) :: nSamples
            character (len=1024) , intent(out):: output_folder
            double precision, dimension (:,:), allocatable, intent(out) :: xMinGlob
            double precision, dimension (:,:), allocatable, intent(out) :: xMaxGlob
            double precision, dimension (:,:), allocatable, intent(out) :: corrL
            integer, dimension(:), allocatable, intent(out) :: corrMod, margiFirst
            integer, dimension(:), allocatable, intent(out) :: seedBase
            character (len=1024), dimension(:), allocatable, intent(out) :: output_name
            double precision, dimension(:), allocatable, intent(out) :: avg, std_dev
            double precision, dimension (:,:), allocatable, intent(out) :: overlap
            integer, dimension(:,:), allocatable, intent(out) :: pointsPerCorrL
            !LOCAL
            integer :: i, j, s
            integer :: fid, fid_2
            integer :: nMat
            double precision, dimension(:,:), allocatable :: bbox_min, bbox_max
            character(len=1024) :: mesh_path, gen_path, absPath
            character(len=1024) :: out_folder
            character(len=1024) :: buffer
            integer, dimension(:), allocatable :: nb_Mat
            integer, dimension(:), allocatable :: nProp_Mat, label_Mat
            integer :: mat, prop_count, mat_Nb, assocMat
            type(property_RF), dimension(:), allocatable :: prop
            integer :: VP_VS_RHO_flag
            double precision :: VP, VS, RHO, VAR
            double precision, dimension(0:2) :: bb_min, bb_max
            logical :: found


            fid_2 = 19

            if(rank == 0) print*, "READING random.spec"
            open (unit = fid_2 , file = trim(adjustL(cur_folder))//"/random.spec", action = 'read', status="old", form="formatted")

                !!1)  Counting number of samples to be generated
                mat = -1
                buffer = getLine(fid_2, '#') !number of materials
                read(buffer,*) nMat
                allocate(nProp_Mat(0:nMat-1))
                allocate(nb_Mat(0:nMat-1))
                do mat = 0, nMat - 1
                    buffer = getLine(fid_2, '#') !Material Number
                    read(buffer,*) nb_Mat(mat)
                    buffer = getLine(fid_2, '#') !VP_VS_RHO and value
                    buffer = getLine(fid_2, '#') !number of properties
                    read(buffer,*) nProp_Mat(mat)
                    do j = 1, nProp_Mat(mat)
                        buffer = getLine(fid_2, '#')
                    end do
                end do

                nSamples = sum(nProp_Mat(:))
                allocate(prop(0:nSamples-1))

                !!2) Reading Samples Properties
                rewind(fid_2)
                prop_count = 0

                buffer = getLine(fid_2, '#') !number of materials
                do mat = 0, nMat - 1
                    buffer = getLine(fid_2, '#') !Material Number
                    buffer = getLine(fid_2, '#') !VP_VS_RHO flag
                    read(buffer,*) VP_VS_RHO_flag
                    if(VP_VS_RHO_flag == 1) then
                        read(buffer,*) VP_VS_RHO_flag, VP, VS, RHO
                        !print*, "VP = ", VP, "VS = ", VS, "RHO = ", RHO
                    end if
                    buffer = getLine(fid_2, '#') !Number of Properties
                    do j = 1, nProp_Mat(mat)
                        buffer = getLine(fid_2, '#')
                        prop(prop_count)%mat = nb_Mat(mat)

                        read(buffer,*) prop(prop_count)%name, &
                                       prop(prop_count)%avg, &
                                       prop(prop_count)%corrMod,   &
                                       prop(prop_count)%corrL(1),  &
                                       prop(prop_count)%corrL(2),  &
                                       prop(prop_count)%corrL(3),  &
                                       prop(prop_count)%margiF,    &
                                       prop(prop_count)%CV,        &
                                       prop(prop_count)%seedStart

                        if(VP_VS_RHO_flag == 1) then
                            if(trim(adjustl(prop(prop_count)%name)) == "Lambda") then
                                prop(prop_count)%avg = 2d0*RHO*VS**2d0
                            else if(trim(adjustl(prop(prop_count)%name)) == "Kappa") then
                                prop(prop_count)%avg = RHO*(VP**2d0 - 4d0*(VS**2d0)/3d0)
                            else if(trim(adjustl(prop(prop_count)%name)) == "Mu") then
                                prop(prop_count)%avg = RHO*VS**2d0
                            else if(trim(adjustl(prop(prop_count)%name)) == "Density") then
                                prop(prop_count)%avg = RHO 
                            end if
                        end if

                        prop_count = prop_count + 1
                    end do
                end do

            close(fid_2)

            if(rank == 0) print*, "READING domains.txt"
            allocate(bbox_min(0:2,0:nMat-1))
            allocate(bbox_max(0:2,0:nMat-1))
            
            open (unit = fid_2 , file = trim(adjustL(cur_folder))//"/domains.txt", action = 'read', status="old", form="formatted")
                
                buffer = getLine(fid_2, '#')
                
                do while (trim(adjustL(buffer)) /= "eof_gl")
                   
                    
                    read(buffer,*) mat_Nb, bb_min(0), bb_min(1), bb_min(2), &
                                          bb_max(0), bb_max(1), bb_max(2)
                    buffer = getLine(fid_2, '#')
                    
                    found = .false.
                    do i = 0, size(nb_Mat)-1
                        if (nb_Mat(i) .eq. mat_Nb) then
                            found = .true.
                            exit
                        end if
                    end do

                    if(.not. found) cycle
                    
                    bbox_min(:,i) = bb_min
                    bbox_max(:,i) = bb_max
                end do

            close(fid_2)

            prop_Count = 0
            do mat = 0, nMat-1
                do j = 0, nProp_Mat(mat)-1
                    mat_Nb = prop(prop_count)%mat
                    prop(prop_Count)%bbox_min = bbox_min(:, mat)
                    prop(prop_Count)%bbox_max = bbox_max(:, mat)
                    prop_Count = prop_count + 1
                end do
            end do

            allocate(output_name(nSamples))
            allocate(xMinGlob(3,nSamples))    
            allocate(xMaxGlob(3,nSamples))
            allocate(corrL(3,nSamples))
            allocate(overlap(3,nSamples))
            allocate(pointsPerCorrL(3,nSamples))
            allocate(corrMod(nSamples))
            allocate(margiFirst(nSamples))
            allocate(avg(nSamples))
            allocate(std_dev(nSamples))
            allocate(seedBase(nSamples))

            output_folder = "./mat/h5"

            do s = 1, nSamples
                output_name(s)               = str_cat("Mat_",num2str(prop(s-1)%mat), &
                                                       "_",prop(s-1)%name)
                xMinGlob(:,s)                = prop(s-1)%bbox_min
                xMaxGlob(:,s)                = prop(s-1)%bbox_max
                corrL(:,s)                   = [prop(s-1)%corrL(1), &
                                                prop(s-1)%corrL(2), &
                                                prop(s-1)%corrL(3)]
                overlap(:,s)                 = [5d0,5d0,5d0]
                pointsPerCorrL(:,s)          = [5,5,5]
                corrMod(s)                   = prop(s-1)%corrMod
                margiFirst(s)                = prop(s-1)%margiF
                avg(s)                       = prop(s-1)%avg
                std_dev(s)                   = prop(s-1)%CV*prop(s-1)%avg
                seedBase(s)                  = prop(s-1)%seedStart
            end do

            if(allocated(prop)) deallocate(prop)
            if(allocated(nProp_Mat)) deallocate(nProp_Mat)
            if(allocated(bbox_max)) deallocate(bbox_max)
            if(allocated(bbox_min)) deallocate(bbox_min)


        end subroutine read_input_SEM
        
    !-------------------------------------------------
    !-------------------------------------------------
    !-------------------------------------------------
    !-------------------------------------------------
    function getLine (fid, comment_Tag) result(nextLine)

        integer,          intent(in) :: fid
        character(len=1), intent(in) :: comment_Tag
        character(len=1024) :: nextLine
        integer :: lineCount = 200, i, stat

        do i = 1, lineCount
            read(fid, fmt="(A)",IOSTAT = stat) nextLine
            nextLine = adjustL(nextLine)
            !print*, "nextLine = ", trim(nextLine)
            if(stat /= 0) then
                nextLine = "eof_gl"
                exit
            else if(nextLine(1:1) /= comment_Tag) then
                !write(*,*) "nextLine = ", trim(nextLine)
                exit
            end if
        end do

    end function getLine
    
    !----------------------------------------------------------------
    !----------------------------------------------------------------
    !----------------------------------------------------------------
    !----------------------------------------------------------------
    function fileExist (filePath) result (fileExists)
        implicit none
        !INPUT
        character(len = *), intent(in) :: filePath
        !OUTPUT
        logical :: fileExists

        inquire(file=filePath, exist=fileExists)
        print*, "file : ", trim(adjustL(filePath))
        print*, "exists?", fileExists

    end function fileExist
end module read_input
