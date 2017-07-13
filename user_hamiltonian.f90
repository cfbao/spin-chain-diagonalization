! Copyright (C) 2012, 2014-2017 - Chenfeng Bao
!
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License; either version 3 of 
! the License, or (at your option) any later version.
! You should have received a copy of the GNU General Public License 
! along with this program; if not, see <http://www.gnu.org/licenses>.

module user_hamiltonian
    use necklaces, only: mod1
    use user_parameters, only: ncell, lcell, nsite, d, nterm, g
    implicit none
contains

pure subroutine h_term_action(alpha_in, s, alpha_out, coeff)
    integer,intent(in) :: alpha_in(:)       ! input state, size = NSITE
    integer,intent(in) :: s                 ! denote this is the s-th term in Hamiltonian
    integer,intent(out) :: alpha_out(:)     ! output state, size = NSITE
    complex(8),intent(out) :: coeff         ! coefficient of output state

! ======================================================================
! ======================= MODIFY BELOW THIS LINE =======================
! ======================================================================
! Example: transverse Ising model
!           H = - sum_i X_i X_i+1 + sum_i Z_i
    integer t
    alpha_out = alpha_in
    if( s <= nsite ) then
! --- implements "- X_i X_i+1" term
        coeff = -1
        t = mod1(s+1,nsite)     ! neighbor index
        alpha_out(s) = mod( alpha_out(s) +1, 2 )
        alpha_out(t) = mod( alpha_out(t)+1, 2 )
    else if( s == nterm ) then
! --- implements "sum Z_i" term
        coeff = 0
        do t = 1, nsite
            coeff = coeff + (-1)**alpha_in(t)
        enddo
    end if
! ======================================================================
! ======================= MODIFY BELOW THIS LINE =======================
! ======================================================================
end subroutine h_term_action
end module user_hamiltonian
