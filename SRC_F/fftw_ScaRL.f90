module fftw_ScaRL
    use constants_ScaRL
    use normal_ScaRL
    use spectra_ScaRL
    use, intrinsic :: iso_c_binding
    !include 'fftw3.f'
    include 'fftw3-mpi.f03'
    !Attention there is a line "include 'fftw3-mpi.f03'" just before this one

contains

subroutine test_FFTw()
    implicit none
    !INPUT
    integer, dimension(3) :: Np
    double precision, dimension(3) :: xRange_in
    integer :: corrMod
    double precision, dimension(3) :: corrL
    integer, dimension(3) :: pointsPerCorrL
    integer ::  seed
    integer :: rank
    !INPUT OUTPUT
    !LOCAL
    integer :: L, M, N
    !type(C_PTR) :: plan
    integer(kind=8) :: plan
    double complex, dimension(:,:,:), allocatable :: Sk_mtx!, gammaK
    !double precision, dimension(:,:,:), allocatable :: r_phase
    double precision :: ampMult
    !double precision, dimension(3) :: delta_k, delta_exp
    double precision, dimension(3) :: xRange
    integer :: ii, jj, kk
    double complex, dimension(:,:,:), allocatable :: c_in
    double complex, dimension(:,:,:), allocatable :: c_out
    !double precision, dimension(3) :: k_max, delta_x, k_cut
    integer, dimension(3) :: sz_c, sz_Sk
    integer :: p_x, p_y, p_z
    double complex :: spec_val
    integer, parameter :: NN=100
    double precision, dimension(NN) :: t, x, w
    double precision :: dt, f
    double complex, dimension(NN) :: in_fft, out_fft
    double complex, dimension(NN) :: Sk_symm_c
    double precision, dimension(NN) :: Sk_half_r, w_half, x_half, t_half
    double precision, dimension(NN) :: in_ifft_r, out_ifft_r
   ! dimension :: in(N), out(N)
    integer*8 planFFT

    dt = 0.01
    f = 1/dt

    do ii = 1, NN
        t(ii) = (ii-1)*dt
        x(ii) = 2*sin(2*PI*10*t(ii)) + 0.5*cos(2*PI*12*t(ii))
        w(ii) = f/dble(NN)*dble(ii-1)
    end do

    do ii = 1, NN/2
        t_half(ii) = 2.*(ii-1)*dt
        w_half(ii) = f/dble(NN/2)*dble(ii-1)
    end do
    
    in_fft = x

   ! call dfftw_plan_dft_1d(planFFT,NN,in_fft,out_fft,FFTW_FORWARD,FFTW_ESTIMATE)
   ! call dfftw_execute_dft(planFFT, in_fft, out_fft)
   ! call dfftw_destroy_plan(planFFT)

    call dfftw_plan_dft_1d(planFFT,NN,in_fft,Sk_symm_c,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(planFFT, in_fft, Sk_symm_c)
    call dfftw_destroy_plan(planFFT)
    
    print *, "time = ["
    do ii = 1, NN
        print *, t(ii)
    end do
    print *, "];"
    
    print *, " x = ["
    do ii = 1, NN
        print *, x(ii)
    end do
    print *, "];"
    
    print *, "w = ["
    do ii = 1, NN
        print *, w(ii)
    end do
    print *, "];"
    
    print *, "t_half = ["
    do ii = 1, NN/2
        print *, t_half(ii)
    end do
    print *, "];"
    
    print *, "w_half = ["
    do ii = 1, NN/2
        print *, w_half(ii)
    end do
    print *, "];"
    
    print *, "Sk_symm_c_R = ["
    do ii = 1, NN
        print *, real(Sk_symm_c(ii))
    end do
    print *, "];"
    
    print *, "Sk_symm_c_I = ["
    do ii = 1, NN
        print *, aimag(Sk_symm_c(ii))
    end do
    print *, "];"

   ! in_fft = out_fft;
   ! call dfftw_plan_dft_1d(planFFT,NN,in_fft,out_fft,FFTW_BACKWARD,FFTW_ESTIMATE)
   ! call dfftw_execute_dft(planFFT, in_fft, out_fft)
   ! call dfftw_destroy_plan(planFFT)
   ! print *, "OUT_iFFT_R = ["
   ! do ii = 1, NN
   !     print *, real(out_fft(ii))
   ! end do
   ! print *, "];"
   ! 
   ! print *, "OUT_iFFT_I = ["
   ! do ii = 1, NN
   !     print *, aimag(out_fft(ii))
   ! end do
   ! print *, "];"
    
    !in_ifft_r = real(Sk_symm_c(1:NN/2))
    Sk_half_r(1:NN/2) = real(Sk_symm_c(1:NN/2))

    
    print *, "Sk_half_r = ["
    do ii = 1, NN/2
        print *, Sk_half_r(ii)
    end do
    print *, "];"
    call dfftw_plan_r2r_1d(planFFT, NN, Sk_half_r, x_half, FFTW_REDFT01, FFTW_ESTIMATE);
    !call dfftw_plan_dft_1d(planFFT,NN/2,in_fft,out_fft,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_r2r(planFFT, Sk_half_r, x_half)
    call dfftw_destroy_plan(planFFT)

    print *, "x_half = ["
    do ii = 1, NN
        print *, x_half(ii)
    end do
    print *, "];"

end subroutine test_FFTw

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
    double complex, dimension(:,:,:), allocatable :: Sk_mtx!, gammaK
    !double precision, dimension(:,:,:), allocatable :: r_phase
    double precision :: ampMult
    !double precision, dimension(3) :: delta_k, delta_exp
    double precision, dimension(3) :: xRange
    integer :: ii, jj, kk
    double complex, dimension(:,:,:), allocatable :: c_in
    double complex, dimension(:,:,:), allocatable :: c_out
    !double precision, dimension(3) :: k_max, delta_x, k_cut
    integer, dimension(3) :: sz_c, sz_Sk
    integer :: p_x, p_y, p_z
    double complex :: spec_val

    integer ( kind = 4 ), parameter :: SS = 10
    real ( kind = 8 ), dimension(2,5) :: r

    xRange = xRange_in/corrL
    sz_c = 2**ceiling(log(dble(2*Np))/log(2d0))
    if(rank ==0)  print*, "Np       = ", Np
    if(rank ==0)  print*, "Size FFT = ", sz_c
    L = sz_c(1)
    M = sz_c(2)
    N = sz_c(3)


    !Build Sk matrix
    call build_Sk_mtx(corrMod, xRange, corrL, Np, sz_c, seed, rank, Sk_mtx)
    sz_Sk = shape(Sk_mtx)/2 

    !Computing data
    data_real_3D = 0d0
    allocate(c_in (sz_c(1),sz_c(2),sz_c(3)))
    
    !c_in(1:sz_Sk(1),1:sz_Sk(2),1:sz_Sk(3)) = Sk_mtx !TESTE +++
    !data_real_3D = real(c_in(1:Np(1),1:Np(2),1:Np(3))) !TESTE

    !if(rank == 0) print*, "maxval(data_real_3D) BEFORE = ", maxval(data_real_3D) 
    !if(rank == 0) print*, "minval(data_real_3D) BEFORE = ", minval(data_real_3D)
    
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
    allocate(c_out(sz_c(1),sz_c(2),sz_c(3)))

    !if(rank == 0) print*, "maxval(c_in) +-- = ", maxval(real(c_in)) 
    !if(rank == 0) print*, "minval(c_in) AFTER = ", minval(real(c_in))
    call dfftw_plan_dft_3d(plan ,L,M,N, c_in,c_out,FFTW_BACKWARD, FFTW_ESTIMATE)
    !plan = fftw_plan_r2r_3d(N,M,L, Sk_mtx, data_real_3d, &
    !       FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, c_in, c_out)
    call dfftw_destroy_plan(plan)
    !if(rank == 0) print*, "maxval(c_out) +-- = ", maxval(real(c_out)) 
    !if(rank == 0) print*, "minval(c_out) AFTER = ", minval(real(c_out))
    data_real_3D = real(c_out(1:Np(1),1:Np(2),1:Np(3)))

    if(rank == 0) print*, "maxval(data_real_3D) AFTER = ", maxval(data_real_3D) 
    if(rank == 0) print*, "minval(data_real_3D) AFTER = ", minval(data_real_3D)

    if(allocated(c_in))   deallocate(c_in)
    if(allocated(c_out))  deallocate(c_out)
    if(allocated(Sk_mtx)) deallocate(Sk_mtx)

end subroutine gen_Std_Gauss_Shino_FFT

!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine gen_std_gauss_FFT(data_real_3D, Np, &
                             xRange, corrL, corrMod, &
                             seed, rank)

    implicit none
    !INPUT
    integer, dimension(3), intent(in) :: Np
    double precision, dimension(3), intent(in) :: xRange
    integer, intent(in) :: corrMod
    double precision, dimension(3), intent(in) :: corrL
    integer, intent(in) ::  seed
    integer, intent(in) :: rank
    !INPUT OUTPUT
    double precision, dimension(:,:,:), intent(inout) :: data_real_3D
    !LOCAL
    integer :: L, M, N
    type(C_PTR) :: plan
    !double precision, dimension(Np(1),Np(2),Np(3)) :: Sk_mtx!, gammaK
    double precision, dimension(:,:,:), allocatable :: Sk_mtx, phiK, gammaK
    double precision, dimension(:,:,:), allocatable :: Sk_mtx_0, data_real_3D_0
    double precision :: ampMult
    !double precision :: k_max
    double precision, dimension(3) :: delta_k, k_max, delta_x
    integer :: ii, jj, kk
    double precision, dimension(:,:,:,:), allocatable :: kPoints_3D

    integer ( kind = 4 ), parameter :: SS = 10
    real ( kind = 8 ), dimension(2,5) :: r
    integer, dimension(:), allocatable :: seedVec
    integer :: seedSz, clock
    double precision :: tmp
    
    allocate(kPoints_3D(3,Np(1)/2,Np(2)/2,Np(3)/2))
    L = Np(1)
    M = Np(2)
    N = Np(3)

    !Build k matrix
    print *, "Build k matrix, proc ", rank
    delta_x(:) = xRange/dble(Np-1)
    k_max(:) = 1.0d0/delta_x(:)
    delta_k(:) = k_max(:)/dble(Np-1)
    !delta_k(:) = 2.0d0*PI/(1.1d0*xRange)
    !k_max(:) = ((Np-1)*delta_k)    
    do kk = 1, Np(3)/2
        do jj = 1, Np(2)/2
            do ii = 1, Np(1)/2
                kPoints_3D(:,ii,jj,kk) = dble([ii,jj,kk]-1)*delta_k(:)/corrL(:)
            end do
        end do
    end do
    
    

    !Build Sk matrix (PSFD)
    print *, "Build Sk matrix, proc ", rank
    allocate(Sk_mtx(Np(1),Np(2),Np(3)))
    select case(corrMod)

    case(cm_GAUSSIAN)
        Sk_mtx(1:Np(1)/2,1:Np(2)/2,1:Np(3)/2) = 1.0d0
        do ii = 1, 3
            Sk_mtx(1:Np(1)/2,1:Np(2)/2,1:Np(3)/2) = Sk_mtx(1:Np(1)/2,1:Np(2)/2,1:Np(3)/2) &
                     * corrL(ii) * &
                      exp(-PI * &
                      (kPoints_3D(ii,:,:,:)**2.0D0))
            !Sk_mtx(1:Np(1)/2,1:Np(2)/2,1:Np(3)/2) = Sk_mtx(1:Np(1)/2,1:Np(2)/2,1:Np(3)/2) &
            !         * corrL(ii) * &
            !          exp(-((kPoints_3D(ii,:,:,:)**2.0D0) * corrL(ii)**2.0d0)/(4.0d0*pi))
        end do
        !Sk_mtx(1+Np(1)/2:,1+Np(2)/2:,1+Np(3)/2:) = 0.d0
      !  Sk_mtx(:,:,:) = 1.0d0/(8d0*PI^3)*&
      !                  exp(-(sum(kPoints_3D(:,:,:,:)**2d0,1)/sum(corrL(:)**2d0))/(4.0d0*pi))
        

  !  case(cm_EXPONENTIAL)
  !      Sk_mtx(:,:,:) = 1.0d0/(8d0*PI^2)*&
  !                      
  !      
  !      do ii = 1, 3
  !          Sk_mtx(:,:,:) = Sk_mtx(:,:,:) * corrL(ii) * &
  !                         exp(-((kPoints_3D(ii,:,:,:)**2.0D0) * corrL(ii)**2.0d0)/(4.0d0*pi))
  !      end do

    end select
    deallocate(kPoints_3D)
    
    !Add random variables    
    print *, "Build random variables, proc ", rank
    allocate(gammaK(Np(1)/2,Np(2)/2,Np(3)/2))
    allocate(phiK(Np(1)/2,Np(2)/2,Np(3)/2))
   ! call r8vec_uniform_01 ( size(phiK), seed+1000, phiK )
   ! call r8vec_normal_01 ( size(gammaK), seed+2000, gammaK )!ERROR TODO find another normal randon number generator
    
    call random_seed(size = seedSz)
    allocate(seedVec(seedSz))
    if(seed >= 0) then
        seedVec = seed + rank
    else
        call system_clock(COUNT=clock)
        seedVec = 100*rank + clock + 37*(/ (ii - 1, ii = 1, seedSz) /) 
    end if 
    call random_seed(put=seedVec)
    call random_number(gammaK) 
    call random_number(phiK) 
    !gammaK  = gammaK -0.5
    gammaK  = 1.0d0 !TEST

    print *, "ADD random variables, proc ", rank
    Sk_mtx(1:Np(1)/2,1:Np(2)/2,1:Np(3)/2)  = gammak*sqrt(Sk_mtx(1:Np(1)/2,1:Np(2)/2,1:Np(3)/2))*cos(2.0D0*PI*phik);
    
    !Process FFT (local)
    !call dfftw_plan_r2r_3d(plan, L, M, N, Sk_mtx, data_real_3d, &
    !       FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
    !call dfftw_execute_r2r(plan, Sk_mtx, data_real_3D)
    !call dfftw_destroy_plan(plan)
    
    plan = fftw_plan_r2r_3d(N,M,L, Sk_mtx, data_real_3d, &
           FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, Sk_mtx, data_real_3D)
    call fftw_destroy_plan(plan)

   ! !TEST
   ! allocate(Sk_mtx_0(0:Np(1)-1,0:Np(2)-1,0:Np(3)-1))
   ! allocate(data_real_3D_0(0:Np(1)-1,0:Np(2)-1,0:Np(3)-1))
   ! do kk = 1, Np(3)
   !     do jj = 1, Np(2)
   !         do ii = 1, Np(1)
   !             Sk_mtx_0(ii-1,jj-1,kk-1) = Sk_mtx(ii,jj,kk)
   !         end do
   !     end do
   ! end do

   ! plan = fftw_plan_r2r_3d(N,M,L, Sk_mtx_0, data_real_3d_0, &
   !        FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
   ! call fftw_execute_r2r(plan, Sk_mtx_0, data_real_3D_0)
   ! call fftw_destroy_plan(plan)
   ! 
   ! do kk = 1, Np(3)
   !     do jj = 1, Np(2)
   !         do ii = 1, Np(1)
   !             data_real_3D(ii,jj,kk) = data_real_3D_0(ii-1,jj-1,kk-1)
   !         end do
   !     end do
   ! end do
   ! if(allocated(Sk_mtx_0)) deallocate(Sk_mtx_0)
   ! if(allocated(data_real_3D_0)) deallocate(data_real_3D_0)
   ! !END_TEST
    
    !data_real_3D = Sk_mtx !TEST
    !data_real_3D = 2.0d0 !TEST

   
    if(allocated(seedVec)) deallocate(seedVec)
    if(allocated(gammaK)) deallocate(gammaK)
    if(allocated(phiK))  deallocate(phiK)

end subroutine gen_Std_Gauss_FFT

!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine gen_Std_Gauss_FFT_MPI(data_real_3D_out, &
                                 xRange, corrL, corrMod, &
                                 seed, rank, L, Np, comm)
    implicit none
    !INPUT
    integer, dimension(3), intent(in) :: L
    integer, intent(in) :: rank, comm
    double precision, dimension(3), intent(in) :: xRange
    integer, intent(in) :: corrMod
    double precision, dimension(3), intent(in) :: corrL
    integer, intent(in) ::  seed
    integer, dimension(3), intent(in) :: Np
    !INPUT OUTPUT
    double precision, dimension(:,:,:), allocatable, intent(out) :: data_real_3D_out
    !LOCAL
    real(C_DOUBLE), pointer :: data_real_3D(:,:,:)
    integer(C_INTPTR_T) :: Np1, Np2, Np3, Np3_offset
    integer(C_INTPTR_T) :: L1,  L2,  L3
    !real(C_DOUBLE), pointer :: data_real_3D(:,:,:)
    type(C_PTR) :: cdata, plan
    integer(C_INTPTR_T) :: alloc_local
    !integer :: sliceSize
    !integer :: kNLocal

    !integer(kind=8) :: i_long, kNCumulated_Init, kNCumulated_End
    integer :: kInit, kEnd
    double precision, dimension(:,:,:,:), allocatable :: kPoints_3D
    double precision, dimension(:,:,:), allocatable :: phiK, gammaK, Sk_mtx
    double precision, dimension(:,:), allocatable :: trash_rand
    double precision :: ampMult
    double precision, dimension(3) :: delta_k, k_max, delta_x
    integer :: ii, jj, kk
    integer, dimension(:), allocatable :: seedVec
    integer :: seedSz, clock
    integer ::nRandVar

    call fftw_mpi_init()

    L1 = L(1); L2 = L(2); L3 = L(3)
    alloc_local = fftw_mpi_local_size_3d(L3, L2, L1, comm, &
                                         Np3, Np3_offset) !FOR MPI
    Np1 = L1; Np2 = L2;

    if(any(Np /= [Np1, Np2, Np3])) stop "Error in Np definition inside FFT_MPI"

    cdata = fftw_alloc_real(alloc_local)
    call c_f_pointer(cdata, data_real_3D, [Np1, Np2, Np3])

    allocate(kPoints_3D(3,Np(1)/2,Np(2)/2,Np(3)))

    !Build k matrix MPI
    print *, "Build k matrix, proc ", rank
    
    kInit = Np3_offset + 1
    kEnd = kInit + Np3 -1
    
    print *, "Proc ", rank, ", kInit = ", kInit
    print *, "Proc ", rank, ", kEnd  = ", kEnd
    print *, "Proc ", rank, ", size(kPoints_3D)  = ", size(kPoints_3D)

    delta_x(:) = xRange/dble(Np-1)
    k_max(:) = 1.0d0/delta_x(:)
    delta_k(:) = k_max(:)/dble(Np-1)
    !delta_k(:) = 2.0d0*PI/(1.1d0*xRange)
    !k_max(:) = ((Np-1)*delta_k)    
    do ii = 1, Np(1)/2
        do jj = 1, Np(2)/2
            do kk = kInit, kEnd
                kPoints_3D(:,ii,jj,kk-kInit+1) = dble([ii,jj,kk]-1)*delta_k(:)/corrL(:)
            end do
        end do
    end do

    !Build Sk matrix MPI
    print *, "Build Sk matrix, proc ", rank
    allocate(Sk_mtx(Np(1),Np(2),Np(3)))
    select case(corrMod)

    case(cm_GAUSSIAN)
        Sk_mtx(1:Np(1)/2,1:Np(2)/2,:) = 1.0d0
        
        do ii = 1, 3
            Sk_mtx(1:Np(1)/2,1:Np(2)/2,:) = Sk_mtx(1:Np(1)/2,1:Np(2)/2,:) * &
                           corrL(ii) * &
                           exp(-PI * &
                           (kPoints_3D(ii,:,:,:)**2.0D0))
        end do

    end select
    
    !allocate(data_real_3D_out(Np1,Np2,Np3)) !FOR DEBUG
    !data_real_3D_out(:,:,:) = Sk_mtx !FOR DEBUG
    !data_real_3D_out(:,:,:) = rank !FOR DEBUG
    !data_real_3D_out(:,:,:) = kPoints_3D(3,:,:,:) !FOR DEBUG
    deallocate(kPoints_3D)



    !Add random variables    
    print *, "Build random variables, proc ", rank
   ! call r8vec_uniform_01 ( size(phiK), seed+1000, phiK )
   ! call r8vec_normal_01 ( size(gammaK), seed+2000, gammaK )!ERROR TODO find another normal randon number generator
    
    call random_seed(size = seedSz)
    allocate(seedVec(seedSz))
    if(seed >= 0) then
        seedVec = seed + rank
        !Putting away the random numbers from others k (that are in others procs)
        allocate(trash_rand(Np(1)/2,Np(2)/2))
        do jj = 1, nRandVar
            do ii = 1, Np3_offset
                call random_number(trash_rand) 
            end do
        end do
        deallocate(trash_rand)
    else
        call system_clock(COUNT=clock)
        seedVec = 100*rank + clock + 37*(/ (ii - 1, ii = 1, seedSz) /) 
    end if 
    call random_seed(put=seedVec)

    allocate(gammaK(Np(1)/2,Np(2)/2,Np(3)))
    allocate(phiK(Np(1)/2,Np(2)/2,Np(3)))
    call random_number(gammaK) 
    call random_number(phiK) 
    !gammaK  = gammaK -0.5
    gammaK  = 1.0 !TEST

    print *, "ADD random variables, proc ", rank
    Sk_mtx(1:Np(1)/2,1:Np(2)/2,:)  = gammak*sqrt(Sk_mtx(1:Np(1)/2,1:Np(2)/2,:))*cos(2.0D0*PI*phik);

    do kk = kInit, kEnd
        if (kk > L(2)/2) then
            Sk_mtx(1+Np(1)/2:,1+Np(2)/2:,kk-kInit+1) = 0.d0
        end if
    end do
    !ampMult = 2.0d0*sqrt(product(MSH%xStep)/((2.0d0)**(dble(RDF%nDim))))

    plan = fftw_mpi_plan_r2r(3, [L3, L2, L1], data_real_3D, data_real_3D, &
                                 comm, [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
    !plan = fftw_mpi_plan_r2r(3, [Np3, Np2, Np1], data_real_3D, data_real_3D, &
    !                             comm, [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
    data_real_3D(:,:,:) = Sk_mtx
    call fftw_mpi_execute_r2r(plan, data_real_3D, data_real_3D)

    write(*,*) "Calculating FFT In rank ", rank
    write(*,*) "shape(data_real_3D) = ", shape(data_real_3D)


    !data_real_3D = data_real_3D*(2.0D0)*sqrt(product(MSH%xStep))
    allocate(data_real_3D_out(Np1,Np2,Np3))
    data_real_3D_out(:,:,:) = data_real_3D
    call fftw_destroy_plan(plan)
    call fftw_free(cdata)

    if(allocated(gammaK)) deallocate(gammaK)
    if(allocated(phik)) deallocate(phik)
    if(allocated(Sk_mtx)) deallocate(Sk_mtx)

end subroutine gen_Std_Gauss_FFT_MPI

end module fftw_ScaRL
