module fftw_ScaRL
    use constants_ScaRL
    use, intrinsic :: iso_c_binding
    !include 'fftw3.f'
    include 'fftw3-mpi.f03'
    !Attention there is a line "include 'fftw3-mpi.f03'" just before this one

contains

!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------


subroutine gen_std_gauss_FFT(data_real_3D, Np, &
                             xRange, corrL, corrMod)

    implicit none
    !INPUT
    integer, dimension(3), intent(in) :: Np
    double precision, dimension(3), intent(in) :: xRange
    integer, intent(in) :: corrMod
    double precision, dimension(3), intent(in) :: corrL
    !INPUT OUTPUT
    double precision, dimension(:,:,:), intent(inout) :: data_real_3D
    !LOCAL
    integer :: L, M, N
    type(C_PTR) :: plan
    double precision, dimension(Np(1),Np(2),Np(3)) :: phiK, gammaK, Sk_mtx
    double precision :: ampMult
    double precision :: k_max
    double precision, dimension(3) :: delta_k
    integer :: ii, jj, kk

    L = Np(1)
    M = Np(2)
    N = Np(3)

    !Build Sk matrix
    select case(corrMod)

    case(cm_GAUSSIAN)
    !k_max = 7.355d0
    delta_k(:) = 2*PI/(xRange+corrL)
    do ii = 1, Np(1)
        do jj = 1, Np(2)
            do kk = 1, Np(3)
                Sk_mtx(ii,jj,kk) = product(corrL) &
                       *exp(-sum((dble([ii,jj,kk]-1)*delta_k(:))**2.0D0 &
                       *(corrL(:)**2.0D0)/(4.0d0*PI)))
            end do
        end do
    end do

    end select
    
    !Add random variables    
    call random_number(gammaK(:,:,:))
    call random_number(phiK(:,:,:))
    gammaK  = gammaK -0.5
    Sk_mtx  = gammak*sqrt(Sk_mtx)*cos(2.0D0*PI*phik);
    
    !Process FFT (local)
    plan = fftw_plan_r2r_3d(N,M,L, Sk_mtx, data_real_3d, &
           FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, Sk_mtx, data_real_3D)
    call fftw_destroy_plan(plan)

end subroutine gen_Std_Gauss_FFT


end module fftw_ScaRL
