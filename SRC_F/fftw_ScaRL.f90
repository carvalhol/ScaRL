module fftw_ScaRL
    use constants_ScaRL
    use normal_ScaRL
    use spectra_ScaRL
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
    double precision :: ampMult
    !double precision :: k_max
    double precision, dimension(3) :: delta_k, k_max
    integer :: ii, jj, kk
    double precision, dimension(:,:,:,:), allocatable :: kPoints_3D

    integer ( kind = 4 ), parameter :: SS = 10
    real ( kind = 8 ), dimension(2,5) :: r
    integer, dimension(:), allocatable :: seedVec
    integer :: seedSz, clock
    double precision :: tmp
    
    allocate(kPoints_3D(3,Np(1),Np(2),Np(3)))
    L = Np(1)
    M = Np(2)
    N = Np(3)

    !Build k matrix
    print *, "Build k matrix, proc ", rank
    delta_k(:) = 2.0d0*PI/(1.1d0*xRange)
    k_max(:) = ((Np-1)*delta_k)    
    do ii = 1, Np(1)
        do jj = 1, Np(2)
            do kk = 1, Np(3)
                kPoints_3D(:,ii,jj,kk) = dble([ii,jj,kk]-1)*delta_k(:)
            end do
        end do
    end do
    
    

    !Build Sk matrix
    print *, "Build Sk matrix, proc ", rank
    allocate(Sk_mtx(Np(1),Np(2),Np(3)))
    select case(corrMod)

    case(cm_GAUSSIAN)
        Sk_mtx(:,:,:) = 1.0d0
        
        do ii = 1, 3
            Sk_mtx(:,:,:) = Sk_mtx(:,:,:) * corrL(ii) * &
                           exp(-((kPoints_3D(ii,:,:,:)**2.0D0) * corrL(ii)**2.0d0)/(4.0d0*pi))
        end do

    end select
    deallocate(kPoints_3D)
    
    !Add random variables    
    print *, "Build random variables, proc ", rank
    allocate(gammaK(Np(1),Np(2),Np(3)))
    allocate(phiK(Np(1),Np(2),Np(3)))
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
    gammaK  = gammaK -0.5

    print *, "ADD random variables, proc ", rank
    Sk_mtx  = gammak*sqrt(Sk_mtx)*cos(2.0D0*PI*phik);
    
    !Process FFT (local)
    plan = fftw_plan_r2r_3d(N,M,L, Sk_mtx, data_real_3d, &
           FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
    call fftw_execute_r2r(plan, Sk_mtx, data_real_3D)
    call fftw_destroy_plan(plan)
   
    if(allocated(seedVec)) deallocate(seedVec)
    if(allocated(gammaK)) deallocate(gammaK)
    if(allocated(phiK))  deallocate(phiK)

end subroutine gen_Std_Gauss_FFT

!subroutine gen_Std_Gauss_FFT_OPT(data_real_3D, Np, &
!                                 xRange, corrL, corrMod, &
!                                 seed)
!    implicit none
!    !INPUT
!    integer, dimension(3), intent(in) :: Np
!    double precision, dimension(3), intent(in) :: xRange
!    integer, intent(in) :: corrMod
!    double precision, dimension(3), intent(in) :: corrL
!    integer, intent(in) ::  seed
!    !INPUT OUTPUT
!    double precision, dimension(:,:,:), intent(inout) :: data_real_3D
!    !LOCAL
!    integer(C_INTPTR_T) :: L, M, N
!    !real(C_DOUBLE), pointer :: data_real_3D(:,:,:)
!    type(C_PTR) :: cdata, plan
!    !integer(C_INTPTR_T) :: alloc_local
!    !integer :: sliceSize
!    !integer :: kNLocal
!
!    double precision, dimension(:), allocatable :: gammaK, phiK
!    integer(kind=8) :: i_long, kNCumulated_Init, kNCumulated_End
!    double precision :: trashNumber
!    integer, dimension(RDF%nDim) :: xNStepGlob
!    double precision :: ampMult
!
!
!    L = Np(1)
!    M = Np(2)
!    N = Np(3)
!
!    call set_kPoints(RDF, MSH%xStep)
!    call set_SkVec(RDF)
!
!
!    kNLocal = size(RDF%kPoints,2)
!
!    allocate(gammaK(kNLocal))
!    allocate(phik(kNLocal))
!
!    gammaK       = gammaK -0.5
!    RDF%SkVec(:) =  gammak*sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);
!    !RDF%SkVec(:) =  sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);
!
!    !TODO Build data_real_3D
!    do ii = 1, Np(1)
!        do jj = 1, Np(2)
!            do kk = 1, Np(3)
!                data_real_3D(ii,jj,kk) = !TODO
!            end do
!        end do
!    end do
!
!    if(allocated(gammaK)) deallocate(gammaK)
!    if(allocated(phik))   deallocate(phik)
!
!    !cdata = fftw_alloc_real(alloc_local)
!    !call c_f_pointer(cdata, data_real_2D, [L,local_M])
!
!    ampMult = 2.0d0*sqrt(product(MSH%xStep)/((2.0d0)**(dble(RDF%nDim))))
!
!    plan = fftw_plan_r2r_3d(N, M, L, data_real_3D, data_real_3D, &
!                             FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01, FFTW_ESTIMATE)
!    !plan = fftw_mpi_plan_r2r(RDF%nDim, [N, M, L], data_real_3D, data_real_3D, &
!    !                         RDF%comm, [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
!    !data_real_3D(:,:,:) = reshape(RDF%SkVec, [L, M, local_LastDim])
!    write(*,*) "Calculating FFT In rang ", RDF%rang
!    !write(*,*) "plan = ", plan
!    write(*,*) "shape(data_real_3D) = ", shape(data_real_3D)
!    !call fftw_mpi_execute_r2r(plan, data_real_3D, data_real_3D)
!    call fftw_execute_r2r(plan, data_real_3D, data_real_3D)
!    write(*,*) "AFTER Calculating FFT In rang ", RDF%rang
!    data_real_3D = data_real_3D*(2.0D0)*sqrt(product(MSH%xStep))
!    !RDF%randField(:,1) = reshape(data_real_3D, [L*M*local_LastDim])
!
!    call fftw_destroy_plan(plan)
!    call fftw_free(cdata)
!
!
!    if(allocated(gammaK)) deallocate(gammaK)
!    if(allocated(phik)) deallocate(phik)
!
!end subroutine gen_Std_Gauss_FFT_OPT

end module fftw_ScaRL
