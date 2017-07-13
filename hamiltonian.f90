module Hamiltonian
    implicit none
    private sortHa
contains

! ============================================================================
! Constructs the Hamiltonian of a translation-invariant spin chain in a 
! K-sector. Periodic boundary condition assumed.
! Result stored in compressed sparse column (CSC) format, lower triangular 
! part only.
! ----------------------------------------------------------------------------
! User must supply H_TERM_ACTION in USER_HAMILTONIAN module to implement a
! custom Hamiltonian.
! ============================================================================
! Arguments
! ----------------------------------------------------------------------------
! HVALUE, HROW, HPNTRB, HPNTRE:
!       (out) {values, row, pointerB, pointerE} of the Hamiltonian in CSC
!             format. See CSC.html for more details.
! K:    (in)  specify in which K-sector the Hamiltonian is to be 
!       constructed
! See CSC.html for details of sparse BLAS CSC matrix storage format
! ============================================================================
subroutine constructH(Hvalue, Hrow, Hpntrb, Hpntre, k)
    use constants
    use user_parameters, only: ncell, nsite, nterm
    use necklaces, only: get_alpha, get_index_x, &
                         norbits, reps, period, kidx, dis2rep
    use user_hamiltonian, only: h_term_action
    ! =========
    ! PARAMETER
    ! =========
    real(8),parameter :: tol = 1.0D-13
    ! =========
    ! ARGUMENTS
    ! =========
    integer,   intent(in)  :: k
    integer,   intent(out) :: Hrow(:), Hpntrb(:), Hpntre(:)
    complex(8),intent(out) :: Hvalue(:)
! ============================================================================
! Local Variables
! ----------------------------------------------------------------------------
! ALPHA_IN:     Initial state, corresponds to a column index
! ALPHA_OUT:    Output state, corresponds to a row index
! S1, S2:       Counter variable for terms in Hamiltonian, range: 1 to NTERM
! COEFFS:       Coefficients in front of output basis vector
! ROWIND:       Row indices of output basis vectors
! DIS:          Number of translations needed for a given spin configuration
!               to become the representative of its equivalent class
! OVERLAP:      Number of ALPHA_OUT's that have the same representative
! HROWT:        Temporary row index
! HCOLUMNT:     Temporary column index
! HVALUET:      Temporary matrix element value
! NHELE:        Counter of matrix elements
! ============================================================================
    integer Hrowt, Hcolumnt, nHele
    complex(8) Hvaluet
    complex(8),allocatable :: coeffs(:)
    integer s1, s2, overlap
    integer alpha_in(nsite), alpha_out(nsite)
    integer,allocatable :: dis(:), rowind(:)
    
    if( size(Hpntrb)/=norbits .or. size(Hpntre)/=norbits .or. size(Hvalue)/=size(Hrow) ) &
        error stop 'Error: inconsistent dimensions in HAMILTONIAN.CONSTRUCTH'
    allocate( coeffs(nterm), dis(nterm), rowind(nterm) )
    ! ==================================
    !          M A I N   L O O P
    ! Loop over different initial states
    !     i.e. different column indices
    ! ==================================
    nHele = 1
    column_loop: do Hcolumnt = 1, norbits
        Hpntrb(Hcolumnt) = nHele
        call get_alpha(reps(Hcolumnt), alpha_in)
        ! ====================================
        ! Loop over all terms in Hamiltonian
        ! constructing states in H*ALPHA_IN
        ! ------------------------------------
        ! H * ket(ALPHA_IN) = H * ket(HCOLUMNT)
        !  = Sum_s1 COEFFS(s1) ket(ALPHA_OUT)
        !  = Sum_s1 COEFFS(s1) ket(ROWIND(s1))
        ! ====================================
        do s1 = 1, nterm
            call h_term_action(alpha_in, s1, alpha_out, coeffs(s1))
            rowind(s1) = kidx(get_index_x(alpha_out))
            dis(s1) = dis2rep(get_index_x(alpha_out))
        enddo
        call sortHa( rowind, dis, coeffs )
        ! =====================================
        ! Loop over different states in H*alpha
        !     i.e. different row indices
        ! =====================================
        s1 = 1
        row_loop: do while( s1 <= nterm )
            Hrowt = rowind(s1)
            ! =============================================
            ! Count # of overlapping H*alpha starting at S1
            ! =============================================
            do overlap = 1, nterm-s1
                if( rowind(s1+overlap)/=rowind(s1) ) exit
            enddo
            !==========================================================
            ! Calculate Hvalue in upper triangle if it might be nonzero
            !==========================================================
            if( mod(k, ncell/period(Hrowt))==0 .and. Hrowt >= Hcolumnt ) then
                ! ============================
                ! Sum over overlapping H*alpha
                ! ============================
                Hvaluet = zero
                do s2 = s1, s1+overlap-1
                    Hvaluet = Hvaluet + coeffs(s2)*exp( onei*2*pi*k*dis(s2)/ncell )
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
    deallocate( coeffs, dis, rowind )
end subroutine constructH

! ==========================================
! Sort ROWIND/DIS/COEFFS in increasing order
!     according to ROWIND
! ==========================================
subroutine sortHa(rowind, dis, coeffs)
    integer,   dimension(:),intent(inout) :: rowind, dis
    complex(8),dimension(:),intent(inout) :: coeffs
    integer n, temp, i, j, key
    complex(8) tempHa
    n = size(rowind)
    if( size(dis)/=n .or. size(coeffs)/=n ) &
        error stop 'Error: inconsistent dimensions in HAMILTONIAN.SORTHA'
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
        
        tempHa = coeffs(i)
        coeffs(i) = coeffs(key)
        coeffs(key) = tempHa
    enddo
end subroutine sortHa

end module Hamiltonian
