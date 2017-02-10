module fftw_ScaRL
    use, intrinsic :: iso_c_binding
    !include 'fftw3.f'
    include 'fftw3-mpi.f03'
    !Attention there is a line "include 'fftw3-mpi.f03'" just before this one

contains

!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------


subroutine gen_std_gauss_FFT(data_real_3D, Np)

    implicit none
    !INPUT
    integer, dimension(3) :: Np
    !INPUT OUTPUT
    double precision, dimension(:,:,:), intent(inout) :: data_real_3D
    !LOCAL
    integer(C_INTPTR_T) :: L, M, N
    integer(C_INTPTR_T) :: local_LastDim
    integer(C_INTPTR_T) :: local_LD_offset
    type(C_PTR) :: cdata, plan
    integer(C_INTPTR_T) :: alloc_local
    integer :: sliceSize
    integer :: kNLocal

    double precision, dimension(Np(1),Np(2),Np(3)) :: phiK, gammaK, k_mtx
    integer(kind=8) :: i_long, kNCumulated_Init, kNCumulated_End
    double precision :: trashNumber
    double precision :: ampMult


    L = Np(1)
    M = Np(2)
    N = Np(3)

    call fftw_mpi_init()

    !alloc_local = fftw_mpi_local_size_3d(N, M, L, RDF%comm, &
    !                                     local_LastDim, local_LD_offset) !FOR MPI
    !cdata = fftw_alloc_real(alloc_local)
    !call c_f_pointer(cdata, data_real_3D, [L, M, local_LastDim])

    !Build k matrix
    k_mtx = 1.0d0
    call random_number(gammaK(:,:,:))
    call random_number(phiK(:,:,:))
    gammaK  = gammaK -0.5

    k_mtx   =  gammak*k_mtx*cos(2.0D0*PI*phik);

    !cdata = fftw_alloc_real(alloc_local)
    !call c_f_pointer(cdata, data_real_2D, [L,local_M])

    !ampMult = 2.0d0*sqrt(product(MSH%xStep)/((2.0d0)**(3d0)))

    plan = fftw_plan_r2r_3d(N, M, L, k_mtx, &
                            FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
   ! plan = fftw_mpi_plan_r2r(RDF%nDim, [N, M, L], data_real_3D, data_real_3D, &
   !                          RDF%comm, [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
    call fftw_mpi_execute_r2r(plan, k_mtx, data_real_3D)

    call fftw_destroy_plan(plan)
    call fftw_free(cdata)

    if(allocated(gammaK)) deallocate(gammaK)
    if(allocated(phik)) deallocate(phik)

end subroutine gen_Std_Gauss_FFT


end module fftw_ScaRL
