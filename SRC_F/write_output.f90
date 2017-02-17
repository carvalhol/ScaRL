module write_output

    use hdf5
    use mpi
    use str_functions
    use hdf5_ScaRL
contains

    !-----------------------------------------------------
    !-----------------------------------------------------
    !-----------------------------------------------------
    !-----------------------------------------------------
    subroutine write_hdf5_multi_proc_3D(xMin, xMax,  &
                                        xMinGlob, xMaxGlob, &
                                        pos_0, pos_N, &
                                        randField_3D, &
                                        HDF5_name,                &
                                        res_folder,           &
                                        rank, nb_procs,           &
                                        comm)

        implicit none

        !INPUTS
        double precision, dimension(3)    , intent(in) :: xMin, xMax
        double precision, dimension(3)    , intent(in) :: xMinGlob, xMaxGlob
        integer         , dimension(3)    , intent(in) :: pos_0, pos_N
        double precision, dimension(:,:,:), intent(in) :: randField_3D
        character(len=*)          :: HDF5_name, res_folder
        integer, intent(in) :: rank, nb_procs, comm
        
        !HDF5 VARIABLES
        character(len=1024)          :: H5_TO_XMF_Path
        integer(HID_T)                 :: file_id       !File identifier
        integer(HID_T)                 :: dset_id       !Dataset identifier
        integer(HID_T)                 :: memspace      ! Dataspace identifier in memory
        integer(HID_T)                 :: filespace
        integer                        :: ds_rank !Dataset rank (number of dimensions)
        integer(HSIZE_T), dimension(3) :: ds_size !Dataset dimensions
        integer                        :: error !Error flag
        integer :: cont = 0

        !LOCAL VARIABLES
        character(len=1024) :: HDF5_path, XMF_path
        double precision, dimension(3) :: xStep
        character(len=1024) :: ds_name
        character(len=1024) :: HDF5_temp
        double precision, dimension(:,:,:), allocatable :: randField_3D_temp
        double precision, dimension(3) :: xMin_temp, xMax_temp
        integer :: i
        integer, dimension(MPI_STATUS_SIZE) :: statut

        if(rank == 0) print*, "HDF5_name = ", trim(HDF5_name)

        if(rank /= 0) then
            call MPI_RECV(cont, 1, MPI_INTEGER, rank-1, 0, comm, statut, error)
            call MPI_SEND(randField_3D, size(randField_3D), MPI_DOUBLE_PRECISION, &
                      0, 0, comm, error)
            call MPI_SEND(xMin, 3, MPI_DOUBLE_PRECISION, &
                      0, 1, comm, error)
            call MPI_SEND(xMax, 3, MPI_DOUBLE_PRECISION, &
                      0, 2, comm, error)
            if(rank /= nb_procs-1) then
            call MPI_SEND(cont, 1, MPI_INTEGER, rank+1, 0, comm, error)
            end if
        else

        !PREPARING ENVIROMENT
        ds_rank = 3
        ds_size = shape(randField_3D)
        ds_name = "sample_1_p"//num2str(rank)
        HDF5_path = str_cat(res_folder,"/",HDF5_name)
        if(rank == 0) print*, "HDF5_path =", trim(HDF5_path)

        !HDF5 WRITING
        call h5open_f(error) ! Initialize FORTRAN interface.
        call h5fcreate_f(HDF5_path, H5F_ACC_TRUNC_F, file_id, error) !NEW file_id
        call h5screate_simple_f(ds_rank, ds_size, filespace, error) !NEW filespace (the size of the whole table)
        call h5dcreate_f(file_id, ds_name, H5T_NATIVE_DOUBLE, filespace, dset_id, error) !NEW dset_id
        call h5screate_simple_f(ds_rank, ds_size, memspace, error)  !NEW memspace

        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                        randField_3D,  &
                        ds_size, error, &
                        file_space_id = filespace, &
                        mem_space_id = memspace) !Write dset, INPUT form = memspace, OUTPUT form = filespace

        call write_h5attr_real_vec(dset_id, "xMin", xMin)
        call write_h5attr_real_vec(dset_id, "xMax", xMax)
        call h5sclose_f(memspace, error) !CLOSE memspace
        call h5dclose_f(dset_id, error) !CLOSE dset_id
        call h5sclose_f(filespace, error) !CLOSE filespace
       
        print*, "nb_procs = ", nb_procs       
        if(nb_procs>1) then
            call MPI_SEND(cont, 1, MPI_INTEGER, rank+1, 0, comm, error)
            allocate(randField_3D_temp(size(randField_3D,1), &
                                   size(randField_3D,2), &
                                   size(randField_3D,3)))
        end if

        do i = 2, nb_procs
        ds_name = "sample_1_p"//num2str(i-1)
        call MPI_RECV(randField_3D_temp, size(randField_3D_temp), MPI_DOUBLE_PRECISION, &
                      i-1, 0, comm, statut, error)
        call MPI_RECV(xMin_temp, 3, MPI_DOUBLE_PRECISION, &
                      i-1, 1, comm, statut, error)
        call MPI_RECV(xMax_temp, 3, MPI_DOUBLE_PRECISION, &
                      i-1, 2, comm, statut, error)
        
        call h5screate_simple_f(ds_rank, ds_size, filespace, error) !NEW filespace (the size of the whole table)
        call h5dcreate_f(file_id, ds_name, H5T_NATIVE_DOUBLE, filespace, dset_id, error) !NEW dset_id
        call h5screate_simple_f(ds_rank, ds_size, memspace, error)  !NEW memspace

        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                        randField_3D_temp,  &
                        ds_size, error, &
                        file_space_id = filespace, &
                        mem_space_id = memspace) !Write dset, INPUT form = memspace, OUTPUT form = filespace

        call write_h5attr_real_vec(dset_id, "xMin", xMin_temp)
        call write_h5attr_real_vec(dset_id, "xMax", xMax_temp)
        call h5sclose_f(memspace, error) !CLOSE memspace
        call h5dclose_f(dset_id, error) !CLOSE dset_id
        call h5sclose_f(filespace, error) !CLOSE filespace
        end do
        
        if(allocated(randField_3D_temp)) deallocate(randField_3D_temp)
        call h5fclose_f(file_id, error) !CLOSE file_id
        call h5close_f(error) ! Close FORTRAN interface

        end if
        

    end subroutine write_hdf5_multi_proc_3D
    
    !-----------------------------------------------------
    !-----------------------------------------------------
    !-----------------------------------------------------
    !-----------------------------------------------------
    subroutine write_hdf5_single_proc_3D(xMin, xMax, &
                                        xMinGlob, xMaxGlob, &
                                        pos_0, pos_N, &
                                        randField_3D, &
                                        HDF5_name, XMF_name,      &
                                        res_folder,           &
                                        rank)

        implicit none

        !INPUTS
        double precision, dimension(3)    , intent(in) :: xMin, xMax
        double precision, dimension(3)    , intent(in) :: xMinGlob, xMaxGlob
        integer         , dimension(3)    , intent(in) :: pos_0, pos_N
        double precision, dimension(:,:,:), intent(in) :: randField_3D
        character(len=*)          :: HDF5_name, XMF_name, res_folder
        integer, intent(in) :: rank
        
        !HDF5 VARIABLES
        character(len=1024)          :: H5_TO_XMF_Path
        integer(HID_T)                 :: file_id       !File identifier
        integer(HID_T)                 :: dset_id       !Dataset identifier
        integer(HID_T)                 :: memspace      ! Dataspace identifier in memory
        integer(HID_T)                 :: filespace
        integer                        :: ds_rank !Dataset rank (number of dimensions)
        integer(HSIZE_T), dimension(3) :: ds_size !Dataset dimensions
        integer                        :: error !Error flag
        
        !LOCAL VARIABLES
        character(len=1024) :: HDF5_path, XMF_path
        double precision, dimension(3) :: xStep
        character(len=1024) :: ds_name

        print*, "HDF5_name = ", trim(HDF5_name)
        print*, "XMF_name = ", trim(XMF_name)

        !PREPARING ENVIROMENT
        ds_rank = 3
        ds_size = shape(randField_3D)
        ds_name = "sample_1_p"//num2str(rank)
        HDF5_path = str_cat(res_folder,"/",HDF5_name)
        XMF_path  = str_cat(res_folder,"/",XMF_name) 

        !HDF5 WRITING
        call h5open_f(error) ! Initialize FORTRAN interface.
        call h5fcreate_f(HDF5_path, H5F_ACC_TRUNC_F, file_id, error) !NEW file_id
        call h5screate_simple_f(ds_rank, ds_size, filespace, error) !NEW filespace (the size of the whole table)
        call h5dcreate_f(file_id, ds_name, H5T_NATIVE_DOUBLE, filespace, dset_id, error) !NEW dset_id
        call h5screate_simple_f(ds_rank, ds_size, memspace, error)  !NEW memspace

        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                        randField_3D,  &
                        ds_size, error, &
                        file_space_id = filespace, &
                        mem_space_id = memspace) !Write dset, INPUT form = memspace, OUTPUT form = filespace

        call write_h5attr_real_vec(dset_id, "xMin", xMin)
        call write_h5attr_real_vec(dset_id, "xMax", xMax)
        call h5sclose_f(memspace, error) !CLOSE memspace
        call h5dclose_f(dset_id, error) !CLOSE dset_id
        call h5sclose_f(filespace, error) !CLOSE filespace
        call h5fclose_f(file_id, error) !CLOSE file_id
        call h5close_f(error) ! Close FORTRAN interface

        !Writing XMF for each part
        if(.false.) then
        call write_XMF_elements(HDF5_name,           &
                                xMin, xMax, shape(randField_3D), &
                                XMF_name, ".", &
                                ".", ds_name)
        end if

    end subroutine write_hdf5_single_proc_3D

    !-----------------------------------------------------------
    !-----------------------------------------------------
    !-----------------------------------------------------
    !-----------------------------------------------------
    !-----------------------------------------------------
    subroutine write_hdf5_multi_proc_3D_1ds(xMin, xMax,  &
                                            xMinGlob, xMaxGlob, &
                                            pos_0, pos_N, &
                                            L, Np, Np_ovlp, &
                                            randField_3D, &
                                            HDF5_name,                &
                                            XMF_name,     &
                                            res_folder,           &
                                            rank, nb_procs,           &
                                            comm)

        implicit none

        !INPUTS
        double precision, dimension(3)    , intent(in) :: xMin, xMax
        double precision, dimension(3)    , intent(in) :: xMinGlob, xMaxGlob
        integer         , dimension(3)    , intent(in) :: pos_0, pos_N
        integer         , dimension(3)    , intent(in) :: L, Np, Np_ovlp
        double precision, dimension(:,:,:), intent(in) :: randField_3D
        character(len=*)          :: HDF5_name, XMF_name, res_folder
        integer, intent(in) :: rank, nb_procs, comm
        
        !HDF5 VARIABLES
        character(len=1024)          :: H5_TO_XMF_Path
        integer(HID_T)                 :: file_id       !File identifier
        integer(HID_T)                 :: dset_id       !Dataset identifier
        integer(HID_T)                 :: memspace      ! Dataspace identifier in memory
        integer(HID_T)                 :: filespace
        integer                        :: ds_rank !Dataset rank (number of dimensions)
        integer(HSIZE_T), dimension(3) :: ds_size !Dataset dimensions
        integer(HSIZE_T), dimension(3) :: local_size !Dataset dimensions
        integer                        :: error !Error flag
        integer(HSSIZE_T), dimension(3) :: offset

        !LOCAL VARIABLES
        character(len=1024) :: HDF5_path, XMF_path
        double precision, dimension(3) :: xStep
        character(len=1024) :: ds_name
        character(len=1024) :: HDF5_temp
        double precision, dimension(:,:,:), allocatable :: randField_3D_temp
        double precision, dimension(3) :: xMin_temp, xMax_temp
        integer :: i
        integer, dimension(MPI_STATUS_SIZE) :: statut
        integer :: cont = 0
        integer, dimension(3) :: pos_0_temp

        if(rank == 0) print*, "write_hdf5_multi_proc_3D_1ds"
        if(rank == 0) print*, "HDF5_name = ", trim(HDF5_name)
        if(any(Np_ovlp < 1)) stop ("When using write_hdf5_multi_proc_3D_1ds overlap should be >= 0")

        if(rank /= 0) then
            call MPI_RECV(cont, 1, MPI_INTEGER, rank-1, 0, comm, statut, error)
            call MPI_SEND(randField_3D, size(randField_3D), MPI_DOUBLE_PRECISION, &
                      0, 0, comm, error)
            call MPI_SEND(xMax, 3, MPI_DOUBLE_PRECISION, &
                      0, 2, comm, error)
            call MPI_SEND(pos_0, 3, MPI_INTEGER, &
                      0, 3, comm, error)
            if(rank /= nb_procs-1) then
            call MPI_SEND(cont, 1, MPI_INTEGER, rank+1, 0, comm, error)
            end if
        else

        !PREPARING ENVIROMENT
        ds_rank = 3
        ds_size = L
        local_size = Np
        ds_name = "sample_1_global"
        HDF5_path = str_cat(res_folder,"/",HDF5_name)
        if(rank == 0) print*, "HDF5_path =", trim(HDF5_path)

        !HDF5 WRITING
        call h5open_f(error) ! Initialize FORTRAN interface.
        call h5fcreate_f(HDF5_path, H5F_ACC_TRUNC_F, file_id, error) !NEW file_id
        call h5screate_simple_f(ds_rank, ds_size, filespace, error) !NEW filespace (the size of the whole table)
        call h5dcreate_f(file_id, ds_name, H5T_NATIVE_DOUBLE, filespace, dset_id, error) !NEW dset_id
        call h5screate_simple_f(ds_rank, local_size, memspace, error)  !NEW memspace
        offset = pos_0 - 1
        call h5sselect_none_f(filespace, error) !DESELECT everything 
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, local_size, error) !SET filespace (to the portion in the hyperslab)

        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                        randField_3D,  &
                        local_size, error, &
                        file_space_id = filespace, &
                        mem_space_id = memspace) !Write dset, INPUT form = memspace, OUTPUT form = filespace

        call write_h5attr_real_vec(dset_id, "xMin", xMin)
      
        if(nb_procs > 1) then 
            call MPI_SEND(cont, 1, MPI_INTEGER, rank+1, 0, comm, error)
            allocate(randField_3D_temp(size(randField_3D,1), &
                                   size(randField_3D,2), &
                                   size(randField_3D,3)))
        end if

        do i = 2, nb_procs
        call MPI_RECV(randField_3D_temp, size(randField_3D_temp), MPI_DOUBLE_PRECISION, &
                      i-1, 0, comm, statut, error)
        call MPI_RECV(xMax_temp, 3, MPI_DOUBLE_PRECISION, &
                      i-1, 2, comm, statut, error)
        call MPI_RECV(pos_0_temp, 3, MPI_DOUBLE_PRECISION, &
                      i-1, 3, comm, statut, error)
        offset = pos_0_temp - 1
        call h5sselect_none_f(filespace, error) !DESELECT everything 
        call h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, offset, local_size, error) !SET filespace (to the portion in the hyperslab)
        
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                        randField_3D_temp,  &
                        local_size, error, &
                        file_space_id = filespace, &
                        mem_space_id = memspace) !Write dset, INPUT form = memspace, OUTPUT form = filespace

      
        if(i == nb_procs) call write_h5attr_real_vec(dset_id, "xMax", xMax_temp)
        end do
        call h5sclose_f(memspace, error) !CLOSE memspace
        call h5dclose_f(dset_id, error) !CLOSE dset_id
        call h5sclose_f(filespace, error) !CLOSE filespace
        
        if(allocated(randField_3D_temp)) deallocate(randField_3D_temp)
        call h5fclose_f(file_id, error) !CLOSE file_id
        call h5close_f(error) ! Close FORTRAN interface

        call write_XMF_elements(HDF5_name,           &
                                xMin, xMax_temp, L, &
                                XMF_name, res_folder, &
                                ".", ds_name)

        end if
        

    end subroutine write_hdf5_multi_proc_3D_1ds
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    subroutine write_XMF_elements(HDF5_name,            & 
                                  xMin, xMax, xNStep,   &
                                  XMF_name, XMF_folder, &
                                  H5_TO_XMF_path, ds_name)

        implicit none

        !INPUTS
        double precision, dimension(3), intent(in) :: xMin, xMax
        integer         , dimension(3), intent(in) :: xNStep
        character(len=*), intent(in) :: HDF5_name
        character(len=*), intent(in) :: XMF_name;
        character(len=*), intent(in) :: XMF_folder
        character(len=*), intent(in) :: H5_TO_XMF_path, ds_name


        !LOCAL VARIABLES
        logical :: DIRECT_OUT = .true.
        integer             :: i, file
        character (len=110) :: dimText;
        double precision, dimension(3) :: xStep

        print*, "-------Writing XMF Mono Processor (Elements)-------"
        print*, "XMF_name   = ", trim(XMF_name)
        print*, "XMF_folder = ", trim(XMF_folder)
        print*, "HDF5_name = ", trim(HDF5_name)
        print*, "H5_TO_XMF_path = ", trim(H5_TO_XMF_path)
        xStep = (xMax-xMin)/(xNStep-1)

        !Writing Number of points in each Dimensions in the reverse order
        dimText = ""
        
        if(DIRECT_OUT) then
            do i = size(xNStep), 1, -1
                dimText = trim(dimText)//" "//trim(num2str(xNStep(i)))
            end do
        else
            do i = 1, size(xNStep)
                dimText = trim(dimText)//" "//trim(num2str(xNStep(i)))
            end do
        end if
        dimText = trim(adjustL(dimText))

        !Building file
        file=21;
        open (unit = file , file = str_cat(XMF_folder,"/",XMF_name), &
              action = 'write')

        write (file,'(A)'      )'<?xml version="1.0" ?>'
        write (file,'(A)'      )'<Xdmf Version="2.0">'
        write (file,'(A)'      )' <Domain>'
        i = 1
        write (file,'(3A)'     )'   <DataItem Name="RF_1" Format="HDF" DataType="Float" Precision="8" &
                                 &Dimensions="',trim(dimText),'">'
        write (file,'(6A)'     )'        ',trim(H5_TO_XMF_path),'/',trim(HDF5_name),':/',trim(adjustL(ds_name))
        write (file,'(A)'      )'   </DataItem>'
        write (file,'(A)'      )'  <Grid GridType="Collection" CollectionType="Spatial">' !Opens the Collection

        write (file,'(A)'     )'   <Grid Name="Group1">'
        write (file,'(3A)'    )'     <Topology TopologyType="3DCoRectMesh" Dimensions="',trim(dimText),'"/>'
        write (file,'(A)'      )'     <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write (file,'(3A)'     )'   <DataItem Name="origin" Format="XML" DataType="Float" &
                                    &Precision="8" Dimensions="',trim(num2str(3)),'">'
        if(DIRECT_OUT) then
             write (file,'(A,F25.10)'      )' ', xMin(3), ' ', xMin(2), ' ', xMin(1)
        else
             write (file,'(A,F25.10)'      )' ', xMin(1), ' ', xMin(2), ' ', xMin(3)
        end if

        write (file,'(A)'      )'   </DataItem>'
        write (file,'(3A)'     )'   <DataItem Name="step" Format="XML" DataType="Float" &
                                    &Precision="8" Dimensions="',trim(num2str(3)),'">'
        if(DIRECT_OUT) then
            write (file,'(A,F25.10)'      )' ', xStep(3), ' ', xStep(2), ' ', xStep(1)
        else
            write (file,'(A,F25.10)'      )' ', xStep(1), ' ', xStep(2), ' ', xStep(3)
        end if

        write (file,'(A)'      )'   </DataItem>'
        write (file,'(A)'     )'     </Geometry>'

        !TODO, DEAL WITH SEVERAL SAMPLES
        do i = 1, 1
            write (file,'(3A)'     )'     <Attribute Name="RF1" Center="Node" AttributeType="Scalar">'
            write (file,'(A)'     )'       <DataItem Reference="XML">'
            write (file,'(3A)'     )'         /Xdmf/Domain/DataItem[@Name="RF_1"]'
            write (file,'(A)'     )'       </DataItem>'
            write (file,'(A)'      )'     </Attribute>'
        end do

        write (file,'(A)'      )'   </Grid>' !END Writing the data of one subdomain
        write (file,'(A)'      )'  </Grid>' !Closes the Collection
        write (file,'(A)'      )' </Domain>'
        write (file,'(A)'      )'</Xdmf>'

        !print*, "------------END Writing XMF Mono Processor (Elements)----------------"

    end subroutine write_XMF_Elements
    
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    subroutine write_XMF_global(globXMF_name, HDF5_list,     & 
                                coord_0_list, xStep, xNStep,  &
                                nb_procs)

        implicit none

        !INPUTS
        character(len=*) :: globXMF_name
        character(len=*), dimension(nb_procs) :: HDF5_list
        double precision, dimension(3, nb_procs) :: coord_0_list
        integer         , dimension(3), intent(in) :: xNStep
        double precision, dimension(3), intent(in) :: xStep
        integer, intent(in) :: nb_procs

        !LOCAL
        character(len=1024) :: XMF_folder="."
        character(len=1024) :: H5_TO_XMF_path=".", ds_name
        logical :: DIRECT_OUT = .true.
        integer             :: i, file
        character (len=110) :: dimText;

        print*, "-------Writing XMF Global-------"

        !Writing Number of points in each Dimensions in the reverse order
        dimText = ""
        
        if(DIRECT_OUT) then
            do i = size(xNStep), 1, -1
                dimText = trim(dimText)//" "//trim(num2str(xNStep(i)))
            end do
        else
            do i = 1, size(xNStep)
                dimText = trim(dimText)//" "//trim(num2str(xNStep(i)))
            end do
        end if
        dimText = trim(adjustL(dimText))

        !Building file
        file=21;
        open (unit = file , file = str_cat(XMF_folder,"/",globXMF_name), &
              action = 'write')

        write (file,'(A)'      )'<?xml version="1.0" ?>'
        write (file,'(A)'      )'<Xdmf Version="2.0">'
        write (file,'(A)'      )' <Domain>'
        
        do i = 1, nb_procs
        write (file,'(3A)'     )'   <DataItem Name="RF_1_p'//trim(num2str(i-1))//'" Format="HDF" DataType="Float" Precision="8" &
                                 &Dimensions="',trim(dimText),'">'
        write (file,'(6A)'     )'        ',trim(H5_TO_XMF_path),'/',trim(HDF5_list(i)),':/sample_1_p'//num2str(i-1)
        write (file,'(A)'      )'   </DataItem>'
        end do

        write (file,'(A)'      )'  <Grid GridType="Collection" CollectionType="Spatial">' !Opens the Collection

        do i = 1, nb_procs
        write (file,'(A)'     )'   <Grid Name="Proc_'//trim(num2str(i-1))//'">'
        write (file,'(3A)'    )'     <Topology TopologyType="3DCoRectMesh" Dimensions="',trim(dimText),'"/>'
        write (file,'(A)'      )'     <Geometry GeometryType="ORIGIN_DXDYDZ">'
        write (file,'(3A)'     )'   <DataItem Name="origin" Format="XML" DataType="Float" &
                                    &Precision="8" Dimensions="',trim(num2str(3)),'">'
        if(DIRECT_OUT) then
             write (file,'(A,F25.10)'      )' ',coord_0_list(3,i), ' ',coord_0_list(2,i), ' ',coord_0_list(1,i)
        else
             write (file,'(A,F25.10)'      )' ',coord_0_list(1,i), ' ',coord_0_list(2,i), ' ',coord_0_list(3,i)
        end if

        write (file,'(A)'      )'   </DataItem>'
        write (file,'(3A)'     )'   <DataItem Name="step" Format="XML" DataType="Float" &
                                    &Precision="8" Dimensions="',trim(num2str(3)),'">'
        if(DIRECT_OUT) then
            write (file,'(A,F25.10)'      )' ', xStep(3), ' ', xStep(2), ' ', xStep(1)
        else
            write (file,'(A,F25.10)'      )' ', xStep(1), ' ', xStep(2), ' ', xStep(3)
        end if

        write (file,'(A)'      )'   </DataItem>'
        write (file,'(A)'     )'     </Geometry>'

        !TODO, DEAL WITH SEVERAL SAMPLES
            write (file,'(3A)'     )'     <Attribute Name="RF_1" Center="Node" AttributeType="Scalar">'
            write (file,'(A)'     )'       <DataItem Reference="XML">'
            write (file,'(3A)'     )'         /Xdmf/Domain/DataItem[@Name="RF_1_p'//trim(num2str(i-1))//'"]'
            write (file,'(A)'     )'       </DataItem>'
            write (file,'(A)'      )'     </Attribute>'

        write (file,'(A)'      )'   </Grid>' !END Writing the data of one subdomain
        end do   
       
        write (file,'(A)'      )'  </Grid>' !Closes the Collection
        write (file,'(A)'      )' </Domain>'
        write (file,'(A)'      )'</Xdmf>'


    end subroutine write_XMF_global
!    !---------------------------------------------------------------
!    !---------------------------------------------------------------
!    !---------------------------------------------------------------
!    !---------------------------------------------------------------
    subroutine mount_hdf5_files(mHDF5_name, nb_procs)
        implicit none
        !INPUT
        character(len=*), intent(in) :: mHDF5_name
        integer, intent(in) :: nb_procs

        !LOCAL
        integer :: f_id, g_id, error
        integer :: ii
        character(len=1024) :: HDF5_name, g_name = "/RF"
        integer :: f_id2

        CALL h5open_f(error)
        
        CALL h5fcreate_f(mHDF5_name, H5F_ACC_TRUNC_F, f_id, error)
        CALL h5gcreate_f(f_id, g_name, g_id, error)
        CALL h5gclose_f(g_id, error)
        CALL h5fclose_f(f_id, error)
        
        CALL h5fopen_f(mHDF5_name, H5F_ACC_RDWR_F, f_id, error)
        
        do ii = 1, nb_procs
            HDF5_name = str_cat("HDF5_proc_", &
                                trim(num2str(ii-1,4)),".h5")
            CALL h5fopen_f(HDF5_name, H5F_ACC_RDWR_F, f_id2, error)
            CALL h5fmount_f (f_id, g_name, f_id2, error)
            CALL h5fclose_f(f_id2, error)
        end do

        CALL h5fclose_f(f_id, error)

        CALL h5close_f(error)


    end subroutine mount_hdf5_files

end module write_output
