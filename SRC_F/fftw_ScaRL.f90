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
                             xRange_in, corrL, &
                             pointsPerCorrL, corrMod, &
                             seed, rank)

    implicit none
    !INPUT
    integer, dimension(3), intent(in) :: Np
    double precision, dimension(3), intent(in) :: xRange_in
    integer, intent(in) :: corrMod
    double precision, dimension(3), intent(in) :: corrL
    integer, dimension(3), intent(in) :: pointsPerCorrL
    integer, intent(in) ::  seed
    integer, intent(in) :: rank
    !INPUT OUTPUT
    double precision, dimension(:,:,:), intent(inout) :: data_real_3D
    !LOCAL
    integer :: L, M, N
    !type(C_PTR) :: plan
    integer(kind=8) :: plan
    double precision, dimension(:,:,:), allocatable :: Sk_mtx!, gammaK
    double precision, dimension(:,:,:), allocatable :: r_phase
    double precision :: ampMult
    double precision, dimension(3) :: delta_k, delta_exp, xRange
    integer :: ii, jj, kk
    double complex, dimension(:,:,:), allocatable :: c_in
    double complex, dimension(:,:,:), allocatable :: c_out
    double precision, dimension(3) :: k_max, delta_x, k_cut
    integer, dimension(3) :: sz_c, sz_Sk
    integer :: p_x, p_y, p_z
    double complex :: spec_val

    integer ( kind = 4 ), parameter :: SS = 10
    real ( kind = 8 ), dimension(2,5) :: r

    xRange = xRange_in/corrL
    sz_c = 2**ceiling(log(dble(2*Np))/log(2d0))
    print*, "Np       = ", Np
    print*, "Size FFT = ", sz_c
    L = sz_c(1)
    M = sz_c(2)
    N = sz_c(3)


    !Build Sk matrix
    select case(corrMod)

    case(cm_GAUSSIAN)
        k_cut = 7.355d0
        delta_x = xRange/dble(Np-1)
        
        k_max = 1/delta_x
        delta_k = k_max/dble(sz_c/2)
        if(any(k_max < k_cut)) then
            print*, "WARNING!! k_max = ", k_max, "< k_cut = ", k_cut
            print*, "Your mesh is too coarse to this representation, truncating k_cut (it may generate an ill-represented field)"
            where(k_max > k_cut) k_cut = k_max
            print*, "New k_cut = ", k_cut
        end if
        sz_Sk = ceiling(k_cut/delta_k)
    
        !Random phase
        allocate(r_phase(sz_Sk(1),sz_Sk(2)*2,sz_Sk(3)*2))
        allocate(Sk_mtx(sz_Sk(1),sz_Sk(2)*2,sz_Sk(3)*2))
        call r8vec_uniform_01 ( size(r_phase), seed+1000000, r_phase)
        r_phase = 2d0*PI*r_phase
        Sk_mtx = cmplx(cos(r_phase),sin(r_phase))
        deallocate(r_phase)
        
        if(rank==0) print*, "delta_x = ", delta_x
        if(rank==0) print*, "k_max   = ", k_max
        if(rank==0) print*, "delta_k = ", delta_k
        if(rank==0) print*, "sz_Sk = ", sz_Sk
        
        delta_exp(:) = ((PI**2d0)/2d0)*((delta_k**2d0))/(corrL**2.0d0)
        if(rank==0) print*, "delta_exp = ", delta_exp
        Sk_mtx(1,:,:)=0d0
        Sk_mtx(:,1,:)=0d0
        Sk_mtx(:,:,1)=0d0
        
        Sk_mtx(:,sz_Sk(2)*2,:)=0d0
        Sk_mtx(:,:,sz_Sk(3)*2)=0d0
        
        do ii = 2, sz_Sk(1)
            do jj = 2, sz_Sk(2)
                do kk = 2, sz_Sk(3)
                    spec_val = sqrt(exp(-(sum((dble([ii,jj,kk]-1)**2d0)*delta_exp(:)))))
                    Sk_mtx(ii,jj,kk) = Sk_mtx(ii,jj,kk)*spec_val !+++
                    Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,kk) = Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,kk)*spec_val !+-+
                    Sk_mtx(ii,jj,(sz_Sk(3)*2)-kk+1) = Sk_mtx(ii,jj,(sz_Sk(3)*2)-kk+1)*spec_val !++-
                    Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,(sz_Sk(3)*2)-kk+1) = Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,(sz_Sk(3)*2)-kk+1)*spec_val !+--
                end do
            end do
        end do
        if(rank==0) print *, "Sk_mtx = ", Sk_mtx(2,2,:)
        Sk_mtx = sqrt(8d0*product(corrL))*Sk_mtx !From the spectra definition
        Sk_mtx = 2d0*sqrt(product(delta_k))*Sk_mtx !Shinozuka formula
    end select
     !Sk_mtx = Sk_mtx/sqrt(dble(product(sz_c))) !FFT normalization
    

    !Computing data
    data_real_3D = 0d0
    allocate(c_in (sz_c(1),sz_c(2),sz_c(3)))
    allocate(c_out(sz_c(1),sz_c(2),sz_c(3)))
    
    !c_in(1:sz_Sk(1),1:sz_Sk(2),1:sz_Sk(3)) = Sk_mtx !TESTE +++
    !data_real_3D = real(c_in(1:Np(1),1:Np(2),1:Np(3))) !TESTE

    if(rank == 0) print*, "maxval(data_real_3D) BEFORE = ", maxval(data_real_3D) 
    if(rank == 0) print*, "minval(data_real_3D) BEFORE = ", minval(data_real_3D)
    
    c_in = 0d0
    !Quadrant(+++)
    c_in(1:sz_Sk(1),&
         1:sz_Sk(2),&
         1:sz_Sk(3)) = Sk_mtx(:,1:sz_Sk(2),1:sz_Sk(3))
    !Quadrant(+-+)
    c_in(1:sz_Sk(1), &
         sz_c(2)-sz_Sk(2)+1:sz_c(2), &
         1:sz_Sk(3)) = Sk_mtx(:,sz_Sk(2)+1:,1:sz_Sk(3))
    !Quadrant(++-)
    c_in(1:sz_Sk(1), &
         1:sz_Sk(2), &
         sz_c(3)-sz_Sk(3)+1:sz_c(3)) = Sk_mtx(:,1:sz_Sk(2), sz_Sk(3)+1:)
    !Quadrant(+--)
    c_in(1:sz_Sk(1), &
         sz_c(2)-sz_Sk(2)+1:sz_c(2), &
         sz_c(3)-sz_Sk(3)+1:sz_c(3)) = Sk_mtx(:,sz_Sk(2)+1:,sz_Sk(3)+1:)
   
    if(allocated(Sk_mtx)) deallocate(Sk_mtx)

    if(rank == 0) print*, "maxval(c_in) +-- = ", maxval(real(c_in)) 
    if(rank == 0) print*, "minval(c_in) AFTER = ", minval(real(c_in))
    call dfftw_plan_dft_3d(plan ,L,M,N, c_in,c_out,FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, c_in, c_out)
    call dfftw_destroy_plan(plan)
    if(rank == 0) print*, "maxval(c_out) +-- = ", maxval(real(c_out)) 
    if(rank == 0) print*, "minval(c_out) AFTER = ", minval(real(c_out))
    data_real_3D = real(c_out(1:Np(1),1:Np(2),1:Np(3)))

    if(rank == 0) print*, "maxval(data_real_3D) +-- = ", maxval(data_real_3D) 
    if(rank == 0) print*, "minval(data_real_3D) AFTER = ", minval(data_real_3D)


   if(allocated(c_in))   deallocate(c_in)
   if(allocated(c_out))  deallocate(c_out)
   if(allocated(Sk_mtx)) deallocate(Sk_mtx)
   if(allocated(phiK))  deallocate(phiK)
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
