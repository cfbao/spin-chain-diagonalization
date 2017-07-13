! Copyright (C) 2012, 2014-2017 - Chenfeng Bao
!
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License; either version 3 of 
! the License, or (at your option) any later version.
! You should have received a copy of the GNU General Public License 
! along with this program; if not, see <http://www.gnu.org/licenses>.

module necklaces
    use user_parameters, only: d, ncell, lcell, nsite
    implicit none
    integer norbits
    integer, allocatable :: reps(:)         ! input index_k output index_x of rep
    integer, allocatable :: period(:)       ! input index_k output period
    integer, allocatable :: kidx(:)         ! input index_x output index_k
    integer, allocatable :: dis2rep(:)      ! input index_x output distance to rep
    private :: get_period, rotate
contains

! ============================================================
! Explicitly find all orbits/necklaces, including their
!   representative, period, k_index, distance_to_rep
!   store in REPS, PERIOD, KIDX, DIS2REP
! This subroutine should be called only once as initialization
! ============================================================
subroutine find_orbits()
    integer alpha(ncell*lcell)
    integer iX, iXt, iK, i
    
    norbits = count_necklaces(ncell, d**lcell)
    allocate( reps(norbits), period(norbits), kidx(d**nsite), dis2rep(d**nsite) )
    kidx = 0
    
    iK = 1
    do iX = 1, d**nsite
        if ( kidx(iX) /= 0 ) cycle
        call get_alpha(iX, alpha)
        reps(iK) = iX
        period(iK) = get_period(alpha)
        ! compute for all states in the orbit
        kidx(iX) = iK
        dis2rep(iX) = 0
        do i = 1, period(iK)-1
            call rotate(alpha)
            iXt = get_index_x(alpha)
            kidx(iXt) = iK
            dis2rep(iXt) = i
        enddo
        iK = iK + 1
    enddo
end subroutine find_orbits

! ======================================
! Calculate the INDEX of given state
! Index is the base-DIM value of a state
! first state is indexed ONE
! ======================================
pure function get_index_x(alpha) result(idx)
    integer,intent(in) :: alpha(:)
    integer idx, i
    idx = 1
    do i = 1, nsite
        idx = idx + alpha(i) * d**(i-1_8)
    enddo
end function get_index_x

! ============================================
! Decipher the given INDEX into explicit state
! return state in ALPHA
! ONE based (i.e. state 0000 has index 1)
! ============================================
pure subroutine get_alpha(idx, alpha)
    integer,intent(in)  :: idx
    integer,intent(out) :: alpha(:)
    integer i, idx0
    idx0 = idx - 1
    do i = 1, nsite
        alpha(i) = mod(idx0, d)
        idx0 = idx0/d
    enddo
end subroutine get_alpha

! ==========================================
! Translate ALPHA by one unit cell
! Def: Translation: alpha(i+lcell) -> alpha(i)
! ==========================================
pure subroutine rotate(alpha)
    integer,intent(inout) :: alpha(:)
    integer i
    integer original(nsite)
    original = alpha
    do i = 1, nsite
        alpha(i) = original( mod1(i+lcell, nsite) )
    enddo
end subroutine

! ========================
! 1-based modulo operation
! result is between [1,n]
! ========================
pure function mod1(m,n)
    integer,intent(in) :: m, n
    integer mod1
    mod1 = modulo(m-1,n)+1
end function mod1

! ========================================================
! return the number of necklaces with n beads and k colors
! requires n >= 1, k >= 2
! ========================================================
pure function count_necklaces(n, k) result(cnt)
    integer,intent(in) :: n, k
    integer cnt
    integer d
    cnt = 0
    do d = 1, n
        if ( mod(n,d) /= 0 ) cycle
        cnt = cnt + phi(d) * k**(n/d)
    enddo
    cnt = cnt / n
    return

contains
    pure function gcd(a, b)
        integer,value :: a, b
        integer gcd, tmp
        do while ( b /= 0 )
            tmp = b
            b = mod(a, b)
            a = tmp
        enddo
        gcd = a
        return
    end function
    
    ! Euler's totient function
    pure function phi(n) result(cnt)
        integer,intent(in) :: n
        integer cnt, i
        cnt = 0
        do i = 1, n
            if ( gcd(i, n) == 1 ) cnt = cnt + 1
        enddo
        return
    end function
end function

! ========================================
! Find the period of the given state ALPHA
! ========================================
pure function get_period(alpha) result(period)
    integer,intent(in) :: alpha(:)
    integer period
    integer i, j
    period = 1
    do
        do i = 1, nsite
            j = mod1( i + lcell*period, nsite )
            if( alpha(i)/=alpha(j) ) exit
        enddo
        if( i > nsite ) return
        period = period + 1
    enddo
end function get_period

end module necklaces
