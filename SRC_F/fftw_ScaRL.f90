module fftw_ScaRL
    use constants_ScaRL
    use normal_ScaRL
    use, intrinsic :: iso_c_binding
    !include 'fftw3.f'
    include 'fftw3-mpi.f03'
    !Attention there is a line "include 'fftw3-mpi.f03'" just before this one

contains
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine gen_std_gauss_Shino_FFT(data_real_3D, Np, &
                             xRange, corrL, corrMod, &
                             seed)

    implicit none
    !INPUT
    integer, dimension(3), intent(in) :: Np
    double precision, dimension(3), intent(in) :: xRange
    integer, intent(in) :: corrMod
    double precision, dimension(3), intent(in) :: corrL
    integer, intent(in) ::  seed
    !INPUT OUTPUT
    double precision, dimension(:,:,:), intent(inout) :: data_real_3D
    !LOCAL
    integer :: L, M, N
    !type(C_PTR) :: plan
    integer(kind=8) :: plan
    double precision, dimension(Np(1),Np(2),Np(3)) :: Sk_mtx!, gammaK
    double precision, dimension(:,:,:,:), allocatable :: r_phase
    double complex, dimension(Np(1),Np(2),Np(3),4) :: phiK
    double precision :: ampMult
    double precision :: k_max
    double precision, dimension(3) :: delta_k, delta_exp
    integer :: ii, jj, kk
    double complex, dimension(:,:,:), allocatable :: c_in
    double complex, dimension(:,:,:), allocatable :: c_out
    integer, dimension(3) :: sz_c
    integer :: p_x, p_y, p_z

    integer ( kind = 4 ), parameter :: SS = 10
    real ( kind = 8 ), dimension(2,5) :: r


    sz_c = 2**ceiling(log(dble(Np))/log(2d0))
    print*, "Size FFT = ", sz_c
    L = sz_c(1)
    M = sz_c(2)
    N = sz_c(3)


    !Build Sk matrix
    Sk_mtx = 0d0
    select case(corrMod)

    case(cm_GAUSSIAN)
        k_max = 7.355d0
        delta_k(:) = 2*PI/xRange
        delta_exp(:) = ((delta_k**2d0) * (corrL**2d0))/(4d0*PI)
        do ii = 2, Np(1)
            do jj = 2, Np(2)
                do kk = 2, Np(3)
                    Sk_mtx(ii,jj,kk) = exp(-sum(dble([ii,jj,kk]-1)*delta_exp(:)))
                    !Sk_mtx(ii,jj,kk) = product(corrL) &
                    !       *exp(-sum((dble([ii,jj,kk]-1)*delta_k(:))**2.0D0 &
                    !       *(corrL(:)**2.0D0)/(4.0d0*PI)))
                end do
            end do
        end do
        Sk_mtx = 2*product(corrL)*product(delta_k)*Sk_mtx
    end select
    !Sk_mtx = Sk_mtx/sqrt(dble(product(sz_c))) !FFT normalization
    !Ps: Sk_mtx supposed symmetric
    !(Sk_mtx(ii,jj,kk) == Sk_mtx(ii,-jj,kk) == Sk_mtx(ii,jj,-kk) == Sk_mtx(ii,-jj,-kk) 
    
    !Random phase
    allocate(r_phase(Np(1),Np(2),Np(3),4))
    call r8vec_uniform_01 ( size(r_phase), seed+1000000, r_phase)
    r_phase = 2d0*PI*r_phase
    phiK = cmplx(cos(r_phase),sin(r_phase))
    deallocate(r_phase)

    !Computing data
    data_real_3D = 0d0
    allocate(c_in (sz_c(1),sz_c(2),sz_c(3)))
    allocate(c_out(sz_c(1),sz_c(2),sz_c(3)))
    
    !Quadrant (+ + +)
    c_in = 0d0
    print*, "phiK(1,1,1,1) = ", phiK(1,1,1,1)
    c_in(1:Np(1),1:Np(2),1:Np(3)) = Sk_mtx*phiK(:,:,:,1)
    call dfftw_plan_dft_3d(plan, L,M,N, c_in,c_out,FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, c_in, c_out)
    call dfftw_destroy_plan(plan)
    data_real_3D = data_real_3D + real(c_out(1:Np(1),1:Np(2),1:Np(3)))




    !Quadrant (+ - +)
    c_in = 0d0
    print*, "phiK(1,1,1,2) = ", phiK(1,1,1,2)
    c_in(1:Np(1), &
         sz_c(2)-Np(2)+1:sz_c(2), &
         1:Np(3)) = Sk_mtx(:,Np(2):1:-1,:)*phiK(:,:,:,2)
    call dfftw_plan_dft_3d(plan, L,M,N, c_in,c_out,FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, c_in, c_out)
    call dfftw_destroy_plan(plan)
    data_real_3D = data_real_3D + real(c_out(1:Np(1),1:Np(2),1:Np(3)))




    !Quadrant (+ + -)
    c_in = 0d0
    print*, "phiK(1,1,1,3) = ", phiK(1,1,1,3)
    c_in(1:Np(1), &
         1:Np(2), &
         sz_c(3)-Np(3)+1:sz_c(3)) = Sk_mtx(:,:,Np(3):1:-1)*phiK(:,:,:,3)
    call dfftw_plan_dft_3d(plan, L,M,N, c_in,c_out,FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, c_in, c_out)
    call dfftw_destroy_plan(plan)
    data_real_3D = data_real_3D + real(c_out(1:Np(1),1:Np(2),1:Np(3)))




    !Quadrant (+ - -)
    c_in = 0d0
    print*, "phiK(1,1,1,4) = ", phiK(1,1,1,4)
    c_in(1:Np(1), &
         sz_c(2)-Np(2)+1:sz_c(2), &
         sz_c(3)-Np(3)+1:sz_c(3)) = Sk_mtx(:,Np(2):1:-1,Np(3):1:-1)*phiK(:,:,:,4)
    call dfftw_plan_dft_3d(plan, L,M,N, c_in,c_out,FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, c_in, c_out)
    call dfftw_destroy_plan(plan)
    data_real_3D = data_real_3D + real(c_out(1:Np(1),1:Np(2),1:Np(3)))


   if(allocated(c_in))  deallocate(c_in)
   if(allocated(c_out)) deallocate(c_out)

end subroutine gen_Std_Gauss_Shino_FFT

!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine gen_std_gauss_FFT(data_real_3D, Np, &
                             xRange, corrL, corrMod, &
                             seed)

    implicit none
    !INPUT
    integer, dimension(3), intent(in) :: Np
    double precision, dimension(3), intent(in) :: xRange
    integer, intent(in) :: corrMod
    double precision, dimension(3), intent(in) :: corrL
    integer, intent(in) ::  seed
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

    integer ( kind = 4 ), parameter :: SS = 10
    real ( kind = 8 ), dimension(2,5) :: r

    L = Np(1)
    M = Np(2)
    N = Np(3)

    !Build Sk matrix
    select case(corrMod)

    case(cm_GAUSSIAN)
    !k_max = 7.355d0
    !delta_k(:) = 2*PI/(xRange+corrL)
    delta_k(:) = 2*PI/(1.1d0*xRange)
    !delta_k(:) = k_max/dble(Np)
    delta_k(:) = ((delta_k**2d0) * (corrL**2d0))/(4d0*PI)
    
    do ii = 1, Np(1)
        do jj = 1, Np(2)
            do kk = 1, Np(3)
                Sk_mtx(ii,jj,kk) = exp(-sum(dble([ii,jj,kk]-1)*delta_k(:)))
                !Sk_mtx(ii,jj,kk) = product(corrL) &
                !       *exp(-sum((dble([ii,jj,kk]-1)*delta_k(:))**2.0D0 &
                !       *(corrL(:)**2.0D0)/(4.0d0*PI)))
            end do
        end do
    end do
    Sk_mtx = product(corrL)*Sk_mtx

    end select
    
    !Add random variables    
    call r8vec_uniform_01 ( size(phiK), seed+1000000, phiK )
    call r8vec_normal_01 ( size(gammaK), seed+2000000, gammaK )
    !call random_number(gammaK(:,:,:))
    !call random_number(phiK(:,:,:))
    !gammaK  = gammaK -0.5
    Sk_mtx  = gammak*sqrt(Sk_mtx)*cos(2.0D0*PI*phik);
    
    !Process FFT (local)
    plan = fftw_plan_r2r_3d(N,M,L, Sk_mtx, data_real_3d, &
           FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, Sk_mtx, data_real_3D)
    call fftw_destroy_plan(plan)

end subroutine gen_Std_Gauss_FFT


end module fftw_ScaRL
