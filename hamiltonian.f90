module Hamiltonian
    implicit none
    private sortHa
contains
! ============================================================================
! This routine constructs the Hamiltonian of a LBP-site-translation-invariant  
! spin-(D-1)/2 chain in a particular K-sector.
! Periodic boundary condition assumed.
! ============================================================================
! Arguments
! ----------------------------------------------------------------------------
! HVALUE:    On output, store numeric value of nonzero matrix elements.
!            Usually the largest array in this routine AND calling program
! HROW:      On output, store row index of nonzero matrix elements
! HCOLUMN:   On output, store column index of nonzero matrix elements
! NHELE:     On output, give number of nonzero matrix elements
! K:         On input, specify in which K-sector the Hamiltonian is to be 
!            constructed
! G:         On input, specify a parameter in the Hamiltonian
!            E.g., a coupling constant.
! ============================================================================
! Local Variables
! ----------------------------------------------------------------------------
! DIS:      Number of translations needed for a given spin configuration to 
!           become the representative of equivalent configurations
! HROWT:    Temporary row index
! HCOLUMNT: Temporary column index
! NTERM:    Number of basis vectors involved after applying the Hamiltonian on
!           a basis vector. Roughly equals the number of nonzero matrix 
!           elements in a column. Usually of the same order as NSITE.
! HVALUET:  Temporary matrix elements' value
! OVERLAP:  Number of ALPHA2's that have the same representative
! ALPHA:    Initial state
! ALPHA2:   Output state
! HALPHA:   Coefficient in front of output basis vector
! ROWIND:   Row indices of all output basis vectors
! S1, S2:   Counter variable for terms of Hamiltonian, range from 1 to NTERM
! ============================================================================
! Specific form of the Hamiltonian must be specified within the implementation
! of the routine. The rule is given by
!         H * ket(ALPHA) = Sum_s1 HVALUE(S1) ket(ALPHA2),   or
!         H * ket(HCOLUMNT) = Sum_s1 HVALUE(S1) ket(ROWIND(S1)).
! More detailed explanations can be found in the note 
!         Constructing Hamiltonians in K-Space
! ============================================================================
subroutine constructH(Hvalue, Hrow, Hpntrb, Hpntre, nterm, k, g)
    use constants
    use necklaces, only: get_alpha, get_index_x, mod1, &
                         nbp, nsite, norbits, reps, period, kidx, dis2rep
    ! =========
    ! PARAMETER
    ! =========
    real(8),parameter :: tol = 1.0D-13
    ! =========
    ! ARGUMENTS
    ! =========
    integer,   intent(in)  :: nterm, k
    integer,   intent(out) :: Hrow(:), Hpntrb(:), Hpntre(:)
    complex(8),intent(out) :: Hvalue(:)
    real(8),   intent(in)  :: g
    ! ===============
    ! LOCAL VARIABLES
    ! ===============
    integer Hrowt, Hcolumnt, nHele
    complex(8) Hvaluet
    complex(8),allocatable :: Halpha(:)
    integer s1, s2, overlap
    integer alpha(nsite), alpha2(nsite)
    integer,allocatable :: dis(:), rowind(:)
    
    if( size(Hpntrb)/=norbits .or. size(Hpntre)/=norbits .or. size(Hvalue)/=size(Hrow) ) &
        error stop 'Error: inconsistent dimensions in HAMILTONIAN.CONSTRUCTH'
    allocate( Halpha(nterm), dis(nterm), rowind(nterm) )
    ! ===========================================
    !            M A I N   L O O P
    ! 1. Construct Hamiltonian in each k-subspace
    ! 2. Count number of elements in each H
    !    store in NHELE
    ! -------------------------------------------
    ! Loop over different initial states
    !     i.e. different column indices
    ! ===========================================
    nHele = 1
    column_loop: do Hcolumnt = 1, norbits
        Hpntrb(Hcolumnt) = nHele
        call get_alpha(reps(Hcolumnt), alpha)
        ! ====================================
        ! Loop over all terms in Hamiltonian
        ! constructing states in H*ALPHA
        ! ------------------------------------
        ! H * ket(ALPHA) = H * ket(HCOLUMNT)
        !  = Sum_s1 HALPHA(S1) ket(ALPHA2)
        !  = Sum_s1 HALPHA(S1) ket(ROWIND(S1))
        ! ====================================
        Halpha = zero
        do s1 = 1, nterm
            alpha2 = alpha
            if( s1 <= nsite ) then
                Halpha(s1) = -1
                alpha2(s1) = mod( alpha2(s1)+1, 2 )
                alpha2(mod1(s1+1,nsite)) = mod( alpha2(mod1(s1+1,nsite))+1, 2 )
            else if( s1 == nterm ) then
                do s2 = 1, nsite
                    Halpha(nterm) = Halpha(nterm) + (-1)**alpha(s2)
                enddo
            end if
            rowind(s1) = kidx(get_index_x(alpha2))
            dis(s1) = dis2rep(get_index_x(alpha2))
        enddo
        call sortHa( rowind, dis, Halpha )
        
        ! =================================
        ! Loop over different states in H*a
        !     i.e. different row indices
        ! =================================
        s1 = 1
        row_loop: do while( s1 <= nterm )
            Hrowt = rowind(s1)
            ! ==========================
            ! Count # of overlapping H*a
            ! ==========================
            do overlap = 1, nterm-s1
                if( rowind(s1+overlap)/=rowind(s1) ) exit
            enddo
            !==========================================================
            ! Calculate Hvalue in upper triangle if it might be nonzero
            !==========================================================
            ! if( mod(k, nBp/period(Hrowt))==0 ) then
            if( mod(k, nBp/period(Hrowt))==0 .and. Hrowt >= Hcolumnt ) then
                ! ========================
                ! Sum over overlapping H*a
                ! ========================
                Hvaluet = zero
                do s2 = s1, s1+overlap-1
                    Hvaluet = Hvaluet + Halpha(s2)*exp( onei*2*pi*k*dis(s2)/nBp )
                enddo
                Hvaluet = Hvaluet * sqrt( 1.0_8*period(Hcolumnt)/period(Hrowt) )
                ! =============================
                ! Save HVALUET if it is nonzero
                ! =============================
                if( abs(Hvaluet) > tol ) then
                    Hvalue (nHele) = Hvaluet
                    Hrow   (nHele) = Hrowt
                    nHele = nHele + 1
                endif
            endif
            ! ======================
            ! Loop counter increment
            ! ======================
            s1 = s1 + overlap
        end do row_loop
        Hpntre(Hcolumnt) = nHele
    end do column_loop
    deallocate( Halpha, dis, rowind )
end subroutine constructH

! ==========================================
! Sort ROWIND/DIS/HALPHA in increasing order
!     according to ROWIND
! ==========================================
subroutine sortHa(rowind, dis, Halpha)
    integer,   dimension(:),intent(inout) :: rowind, dis
    complex(8),dimension(:),intent(inout) :: Halpha
    integer n, temp, i, j, key
    complex(8) tempHa
    n = size(rowind)
    if( size(dis)/=n .or. size(Halpha)/=n ) error stop 'Error: inconsistent dimensions in HAMILTONIAN.SORTHA'
    ! =========================================
    ! Selection sort
    ! I: size of selected array after the loop
    ! J: pointer in unsorted array for scanning
    ! =========================================
    do i = 1, n-1
        key = i
        do j = i+1, n
            if( rowind(j)<rowind(key) ) key = j
        enddo
        
        temp = rowind(i)
        rowind(i) = rowind(key)
        rowind(key) = temp
        
        temp = dis(i)
        dis(i) = dis(key)
        dis(key) = temp
        
        tempHa = Halpha(i)
        Halpha(i) = Halpha(key)
        Halpha(key) = tempHa
    enddo
end subroutine sortHa

end module Hamiltonian
