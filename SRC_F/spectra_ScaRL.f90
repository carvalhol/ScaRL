module spectra_ScaRL
    use constants_ScaRL
    use normal_ScaRL

contains
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
!-------------------------------------------------------------
subroutine build_Sk_mtx(corrMod, xRange, corrL, Np, sz_c, seed, rank, Sk_mtx)

    implicit none
    !INPUT
    integer, dimension(3), intent(in) :: Np
    double precision, dimension(3), intent(in) :: xRange, corrL
    integer, intent(in) :: corrMod
    integer, intent(in) ::  seed
    integer, intent(in) :: rank
    !INPUT OUTPUT
    double complex, dimension(:,:,:), allocatable :: Sk_mtx!, gammaK
    !LOCAL
    double precision, dimension(:,:,:), allocatable :: r_phase
    double precision :: ampMult
    double precision, dimension(3) :: delta_k, delta_exp
    integer :: ii, jj, kk
    double precision, dimension(3) :: k_max, delta_x, k_cut
    integer, dimension(3) :: sz_c, sz_Sk
    integer :: p_x, p_y, p_z
    double complex :: spec_val

    integer ( kind = 4 ), parameter :: SS = 10
    real ( kind = 8 ), dimension(2,5) :: r

    double precision,parameter :: nu=0.1d0 ! < 1.0!!
    double precision :: vm
    double precision, dimension(1) :: bes_i,dbes_i,bes_k,dbes_k
    double precision, dimension(:), allocatable :: kk2
    
    if(rank ==0)  print*, "Building Sk_mtx"

    !Defining size    
    select case(corrMod)
    
    case(cm_GAUSSIAN)
        k_cut = 7.355d0
        delta_x = xRange/dble(Np-1)
        
        k_max = 1/delta_x
        delta_k = k_max/dble(sz_c/2)
        if(any(k_max < k_cut) .and. rank == 0) then
            print*, "WARNING!! k_max = ", k_max, "< k_cut = ", k_cut
            print*, "Your mesh is too coarse to this representation, truncating k_cut (it may generate an ill-represented field)"
            where(k_max > k_cut) k_cut = k_max
            print*, "New k_cut = ", k_cut
        end if
        sz_Sk = ceiling(k_cut/delta_k)

    case(cm_VON_KARMAN)
        k_cut = 7.355d0
        delta_x = xRange/dble(Np-1)
        
        k_max = 1/delta_x
        delta_k = k_max/dble(sz_c/2)
        if(any(k_max < k_cut) .and. rank == 0) then
            print*, "WARNING!! k_max = ", k_max, "< k_cut = ", k_cut
            print*, "Your mesh is too coarse to this representation, truncating k_cut (it may generate an ill-represented field)"
            where(k_max > k_cut) k_cut = k_max
            print*, "New k_cut = ", k_cut
        end if
        sz_Sk = ceiling(k_cut/delta_k)

    end select

    allocate(r_phase(sz_Sk(1),sz_Sk(2)*2,sz_Sk(3)*2))
    allocate(Sk_mtx(sz_Sk(1),sz_Sk(2)*2,sz_Sk(3)*2))
    
    !Random phase
    call r8vec_uniform_01 ( size(r_phase), seed+1000000, r_phase)
    r_phase = 2d0*PI*r_phase
    Sk_mtx = cmplx(cos(r_phase),sin(r_phase))
    Sk_mtx(1,:,:)=0d0
    Sk_mtx(:,1,:)=0d0
    Sk_mtx(:,:,1)=0d0
    Sk_mtx(:,sz_Sk(2)*2,:)=0d0
    Sk_mtx(:,:,sz_Sk(3)*2)=0d0
    deallocate(r_phase)

    if(rank==0) print*, "delta_x = ", delta_x
    if(rank==0) print*, "k_max   = ", k_max
    if(rank==0) print*, "delta_k = ", delta_k
    if(rank==0) print*, "sz_Sk = ", sz_Sk

    
    !Calculate Matrix
    select case(corrMod)
    
    case(cm_GAUSSIAN)
        !delta_exp(:) = ((PI**2d0)/2d0)*((delta_k**2d0))/(corrL**2.0d0)
        delta_exp(:) = ((PI**2d0)/2d0)*((delta_k**2d0))
        if(rank==0) print*, "delta_exp = ", delta_exp
        
        do ii = 2, sz_Sk(1)
            do jj = 2, sz_Sk(2)
                do kk = 2, sz_Sk(3)
                    !spec_val = sqrt(exp(-(sum(((dble([ii,jj,kk]-1)**2d0)*delta_exp(:))/corrL**2d0))))
                    spec_val = sqrt(exp(-(sum(((dble([ii,jj,kk]-1)**2d0)*delta_exp(:))))))
                    Sk_mtx(ii,jj,kk) = Sk_mtx(ii,jj,kk)*spec_val !+++
                    Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,kk) = Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,kk)*spec_val !+-+
                    Sk_mtx(ii,jj,(sz_Sk(3)*2)-kk+1) = Sk_mtx(ii,jj,(sz_Sk(3)*2)-kk+1)*spec_val !++-
                    Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,(sz_Sk(3)*2)-kk+1) = Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,(sz_Sk(3)*2)-kk+1)*spec_val !+--
                end do
            end do
        end do
        !if(rank==0) print *, "Sk_mtx = ", Sk_mtx(2,2,:)
        !Sk_mtx = sqrt(8d0*product(corrL))*Sk_mtx !From the spectra definition
        Sk_mtx = sqrt(8d0)*Sk_mtx !From the spectra definition
        Sk_mtx = 2d0*sqrt(product(delta_k))*Sk_mtx !Shinozuka formula
    
    case(cm_VON_KARMAN)
       ! ! BESSEL FUNCTION
       ! call ikv(nu,1.0d-10,vm,bes_i,dbes_i,bes_k,dbes_k)
       ! ampMult = (corrL(1)**2+corrL(2)**2+corrL(3)**2)*4.0d0*pi*nu/(bes_k(1))
       ! delta_exp(:) = (delta_k*corrL)**2d0
       ! if(rank==0) print*, "delta_exp = ", delta_exp
       ! 
       ! do ii = 2, sz_Sk(1)
       !     do jj = 2, sz_Sk(2)
       !         do kk = 2, sz_Sk(3)
       !             spec_val = ampMult*(1.0d0+(delta_exp(1)*dble(ii**2))+(delta_exp(2)*dble(jj**2))+(delta_exp(3)*dble(kk**2)))**(nu+1.5d0)

       !             Sk_mtx(ii,jj,kk) = Sk_mtx(ii,jj,kk)*spec_val !+++
       !             Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,kk) = Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,kk)*spec_val !+-+
       !             Sk_mtx(ii,jj,(sz_Sk(3)*2)-kk+1) = Sk_mtx(ii,jj,(sz_Sk(3)*2)-kk+1)*spec_val !++-
       !             Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,(sz_Sk(3)*2)-kk+1) = Sk_mtx(ii,(sz_Sk(2)*2)-jj+1,(sz_Sk(3)*2)-kk+1)*spec_val !+--
       !         end do
       !     end do
       ! end do
    
    end select


end subroutine build_Sk_mtx


end module spectra_ScaRL
