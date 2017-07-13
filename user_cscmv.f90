! Copyright (C) 2012, 2014-2017 - Chenfeng Bao
!
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License; either version 3 of 
! the License, or (at your option) any later version.
! You should have received a copy of the GNU General Public License 
! along with this program; if not, see <http://www.gnu.org/licenses>.

module user_cscmv
    use constants, only: one, zero
    implicit none
contains

! ========================================================================
! Wrapper for a custom matrix-vector multiplication routine for a complex
! Hermitian matrix in CSC format, where only the lower triangular part is
! stored.
! ------------------------------------------------------------------------
! Arguments
! N:    (in)  dimension of the matrix
! VALUES, ROWS, POINTERB, POINTERE:
!       (in)  arrays specifying a sparse matrix 'OP' in CSC format
!             see CSC.html for more details
! X:    (in)  input vector
! Y:    (out) output vector, should equal OP * X on return
! ------------------------------------------------------------------------
! This default implementation uses a MKL sparse BLAS level 2 routine.
! Modify this file only if MKL is not a viable option for you.
! ========================================================================
subroutine cscmv(n, values, rows, pointerB, pointerE, x, y)
    integer, intent(in) :: n        ! dimension of matrix
    integer, intent(in),dimension(:) :: rows, pointerB, pointerE
    complex(8),intent(in),dimension(:) :: values
    complex(8),intent(in), dimension(*) :: x
    complex(8),intent(out),dimension(*) :: y
! ======================================================================
! ======================= MODIFY BELOW THIS LINE =======================
! ======================================================================
    character :: matdescra(6) = ['H', 'L', 'N', 'F', 'N', 'N']
    call mkl_zcscmv('N', n, n, one, matdescra, &
                    values, rows, pointerB, pointerE, x, zero, y)
! ======================================================================
! ======================= MODIFY ABOVE THIS LINE =======================
! ======================================================================
end subroutine cscmv

end module user_cscmv
