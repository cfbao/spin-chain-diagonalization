module cycleUtility
    implicit none
    integer d, nbp, lbp, nsite
    integer,    allocatable :: period(:)
    integer(8), allocatable :: config(:)
    private :: is_alpha_min, get_period, update_alpha, get_index_x
contains

! =====================================================
! 1. Find indices/periods of representative of cycles
!    store them in array CONFIGS/PERIODS
! 2. Adjust CONFIGS/PERIODS to exact size (NCONFIG)
! =====================================================
subroutine find_cycles(nbp_in, lbp_in, d_in)
    integer,intent(in)  :: nbp_in, lbp_in, d_in
    integer nconfig
    integer flag, alpha(nbp_in*lbp_in)
    integer,   allocatable :: ptmp(:)
    integer(8),allocatable :: ctmp(:)
    
    nbp = nbp_in
    lbp = lbp_in
    nsite = nbp * lbp
    d = d_in
    nconfig = int((1.0_8*d)**nsite / nbp + (1.0_8*d)**(nsite/2) / nbp * 2 + 2)
    
    if( allocated(period) .and. size(period)<nconfig ) deallocate(period)
    if( allocated(config) .and. size(config)<nconfig ) deallocate(config)
    if( .not.allocated(period) ) allocate( period(nconfig) )
    if( .not.allocated(config) ) allocate( config(nconfig) )
    alpha = 0
    nconfig = 0
    do
        period(nconfig+1) = get_period(alpha)
        if( is_alpha_min(alpha, period(nconfig+1)) ) then
            nconfig = nconfig + 1
            config(nconfig) = get_index_x(alpha)
        endif
        call update_alpha(alpha, flag)
        if(flag == -1) exit
    end do
    if( size(period)/=nconfig ) then
        allocate( ptmp(nconfig) )
        ptmp = period(1:nconfig)
        call move_alloc(ptmp, period)
    endif
    if( size(config)/=nconfig ) then
        allocate( ctmp(nconfig) )
        ctmp = config(1:nconfig)
        call move_alloc(ctmp, config)
    endif
end subroutine find_cycles

! ======================================
! Calculate the INDEX of given state
! Index is the base-DIM value of a state
! first state is indexed ZERO
! ======================================
pure function get_index_x(alpha) result(idx)
    integer,intent(in) :: alpha(:)
    integer(8) idx
    integer nsite, i
    nsite = size(alpha)
    idx = 0
    do i = 1, nsite
        idx = idx + alpha(i) * d**(i-1_8)
    enddo
end function get_index_x

! ==================================================
! Find the K-index of the state (ALPHA)
! Return the value in KINDEX
! K-index is the index in K-subspace (1 to NCONFIG)
! The given state (ALPHA) must be the rep of a cycle
! ==================================================
pure function get_index_k(alpha) result(Kindex)
    integer,   intent(in) :: alpha(:)
    integer :: Kindex
    integer nsite, nconfig, head, tail, mid
    integer(8) idx
    nsite = size(alpha)
    nconfig = size(config)
    
    idx = get_index_x(alpha)
    head = 1
    tail = nconfig
    mid  = head + (tail-head) / 2
    do
        if( idx == config(mid) ) then
            Kindex = mid
            return
        else if( idx > config(mid) ) then
            if( idx == config(tail) ) then
                Kindex = tail
                return
            endif
            head = mid
        else
            if( idx == config(head) ) then
                Kindex = head
                return
            endif
            tail = mid
        endif
        mid = head + (tail-head) / 2
    enddo
end function get_index_k

! ============================================
! Decipher the given INDEX into explicit state
! return state in ALPHA
! ZERO based (i.e. state 0000 has index 0)
! ============================================
pure subroutine get_alpha(idx, alpha)
    integer(8),intent(in)  :: idx
    integer,   intent(out) :: alpha(:)
    integer nsite, i
    integer(8) dummy
    nsite = size(alpha)
    dummy = idx
    do i = 1, nsite
        alpha(i) = mod(dummy, int(d,8))
        dummy = dummy/d
    enddo
end subroutine get_alpha

! =============================================
! 1. Find the representative of the cycle which
!        the given state (ALPHA) belongs to
!    Return the representative in ALPHA
! 2. Find the number of translations needed to
!        become the representative
!    Return this number in DIS
! Def: Translation: alpha(i+LBP) -> alpha(i)
! =============================================
pure subroutine get_rep(alpha, dis)
    integer,intent(inout) :: alpha(:)
    integer,intent (out)  :: dis
    integer nsite, msp, rot, i, j1, j2
    integer,allocatable :: original(:)
    nsite = size(alpha)
    msp = 0
    do rot = lbp, nsite-lbp, lbp
        do i = nsite, 1, -1
            j1 = mod1( i+msp, nsite )
            j2 = mod1( i+rot, nsite )
            if( alpha(j2) < alpha(j1) ) then
                msp = rot
                exit
            else if( alpha(j1) < alpha(j2) ) then
                exit
            endif
        enddo
    enddo
    dis = msp / lbp
    allocate( original(nsite) )
    original = alpha
    do i = 1, nsite
        j1 = mod1( msp+i, nsite )
        alpha(i) = original(j1)
    enddo
    deallocate( original )
end subroutine get_rep

! ========================
! 1-based modulo operation
! result is between [1,n]
! ========================
pure function mod1(m,n)
    integer,intent(in) :: m, n
    integer mod1
    mod1 = modulo(m-1,n)+1
end function mod1

! =========================================================
! Check if the given alpha is the representative of a cycle
! If yes, return .TRUE.; if not, return .FALSE.
! A representative of a cycle is the state (alpha) with the
!     smallest index in the cycle
! =========================================================
pure function is_alpha_min(alpha, period) result(isMin)
    integer,intent(in) :: alpha(:), period
    logical isMin
    integer rot, nsite, i, j
    nsite = size(alpha)
    do rot = 1, period-1
        do i = nsite, 1, -1
            j = mod1( i + lbp*rot, nsite )
            if( alpha(j) < alpha(i) ) then
                isMin = .false.
                return
            else if( alpha(i) < alpha(j) ) then
                exit
            endif
        enddo
    enddo
    isMin = .true.
    return
end function is_alpha_min

! ========================================
! Find the period of the given state ALPHA
! ========================================
pure function get_period(alpha) result(period)
    integer,intent(in) :: alpha(:)
    integer period
    integer nsite, i, j
    nsite = size(alpha)
    period = 1
    do
        do i = 1, nsite
            j = mod1( i + lbp*period, nsite )
            if( alpha(i)/=alpha(j) ) exit
        enddo
        if( i > nsite ) return
        period = period + 1
    enddo
end function get_period

! =================================================
! Update ALPHA similar to binary increment
! Return the updated alpha in the original alpha
! ALPHA(1) represents the smallest digit
! If ALPHA is already the largest, return FLAG = -1
! =================================================
pure subroutine update_alpha(alpha, flag)
    integer,intent(inout) :: alpha(:)
    integer,intent (out)  :: flag
    integer i, j, nsite
    nsite = size(alpha)
    do i = 1, nsite
        if(alpha(i) < d-1 ) exit
    enddo
    if( i > nsite ) then
        flag = -1
        return
    else
        flag = 1
    endif
    alpha(i) = alpha(i) + 1
    do j = 1, i-1
        alpha(j) = 0
    enddo
end subroutine update_alpha

end module cycleUtility
