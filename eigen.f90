! Copyright (C) 2012, 2014-2017 - Chenfeng Bao
!
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License; either version 3 of 
! the License, or (at your option) any later version.
! You should have received a copy of the GNU General Public License 
! along with this program; if not, see <http://www.gnu.org/licenses>.

module eigen
    implicit none
contains

! ==========================================================================
! This routine is an interface to ARPACK. It computes selected eigenvalues 
! (and corresponding eigenvectors if requested) of a complex Hermitian 
! matrix in compressed sparse column (CSC) format, where only the lower 
! triangular part is stored.
! It is modified from the ARPACK example program ZNDRV1.F.
! --------------------------------------------------------------------------
! User may wish to supply a custom subroutine for CSC matrix-vector 
! multiplication if using MKL is not a viable option.
! ==========================================================================
! Arguments
! --------------------------------------------------------------------------
! HVALUE, HROW, HPNTRB, HPNTRE:
!           (in)  {values, row, pointerB, pointerE} of the matrix in CSC format
!                 see CSC.html for more details.
! NEV:      (in)  number of eigenvalues to be computed
! D:        (out) computed eigenvalues
! V:        (out) workspace and optionally store computed eigenvectors on output
! RVEC:     (in)  specify whether eigenvectors are to be computed
! FILCONV:  (in)  unit specifier for log file that stores convergence info
! FILRES:   (in)  unit specifier for log file that stores eigenvalues and 
!                 (if RVEC = .true.) residual info
! ==========================================================================
subroutine eig(Hvalue,Hrow,Hpntrb,Hpntre,nev, D,V, rvec, filconv, filres)
    use user_cscmv, only: cscmv
!   %--------------%
!   |  Arguments   |
!   %--------------%
    integer Hrow(:), Hpntrb(:), Hpntre(:), nev, filconv, filres
    logical rvec
    complex(8) Hvalue(:), D(:), V(:,:)
    
!   %---------------%
!   | Local Scalars |
!   %---------------%
    character bmat*1, which*2
    integer n, ncv, ldv, ido, lworkl, info, j, ierr, nconv, maxitr, ishfts, mode, filerr
    complex(8) sigma
    real(8) tol
    
!   %--------------%
!   | Local Arrays |
!   %--------------%
    integer iparam(11), ipntr(14)
    logical,allocatable :: select(:)
    complex(8),allocatable,dimension(:) :: ax, workd, workev, resid, workl
    real(8),allocatable :: rwork(:), rd(:,:)
    
!   %-----------------------------%
!   | BLAS & LAPACK routines used |
!   %-----------------------------%
    real(8) dznrm2 , dlapy2
    external dznrm2 , zaxpy , dlapy2
    
!   %-----------------------%
!   | Executable Statements |
!   %-----------------------%
    n = size(Hpntrb)
    ncv = size(D)
    ldv = n
    allocate( select(ncv) )
    allocate( ax(n), workd(3*n), workev(3*ncv), resid(n), workl(3*ncv*ncv+5*ncv) )
    allocate( rwork(ncv), rd(ncv,3) )
    
!   %--------------------------------------------------%
!   | The number N is the dimension of the             |
!   | matrix.  A standard eigenvalue problem is        |
!   | solved (BMAT = 'I').  NEV is the number of       |
!   | eigenvalues to be approximated. The user can     |
!   | modify N, NEV, NCV, WHICH to solve problems of   |
!   | different sizes, and to get different parts of   |
!   | the spectrum.  However, The following            |
!   | conditions must be satisfied:                    |
!   |           NEV + 2 <= NCV <= N                    | 
!   %--------------------------------------------------% 
    bmat  = 'I'
    which = 'SR'
    
!   %---------------------------------------------------%
!   | The work array WORKL is used in ZNAUPD  as        | 
!   | workspace.  Its dimension LWORKL is set as        |
!   | illustrated below.  The parameter TOL determines  |
!   | the stopping criterion. If TOL<=0, machine        |
!   | precision is used.  The variable IDO is used for  |
!   | reverse communication, and is initially set to 0. |
!   | Setting INFO=0 indicates that a random vector is  |
!   | generated to start the ARNOLDI iteration.         | 
!   %---------------------------------------------------%
    lworkl = 3*ncv**2 + 5*ncv
    tol = 0.0_8
    ido = 0
    info = 0
    
!   %---------------------------------------------------%
!   | This program uses exact shift with respect to     |
!   | the current Hessenberg matrix (IPARAM(1) = 1).    |
!   | IPARAM(3) specifies the maximum number of Arnoldi |
!   | iterations allowed.  Mode 1 of ZNAUPD  is used    |
!   | (IPARAM(7) = 1). All these options can be changed |
!   | by the user. For details see the documentation in |
!   | ZNAUPD .                                          |
!   %---------------------------------------------------%
    ishfts = 1
    maxitr = 300
    mode = 1
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode
    
!   %-------------------------------------------%
!   | M A I N   L O O P (Reverse communication) | 
!   %-------------------------------------------%
    do
!       %---------------------------------------------%
!       | Repeatedly call the routine ZNAUPD  and take|
!       | actions indicated by parameter IDO until    |
!       | either convergence is indicated or MAXITR   |
!       | has been exceeded.                          |
!       %---------------------------------------------%
        call znaupd(ido, bmat, n, which, nev, tol, resid, ncv, v, ldv, &
                    iparam, ipntr, workd, workl, lworkl, rwork, info )
        if( ido /= -1 .and. ido /= 1 ) exit
!       %-------------------------------------------%
!       | Perform matrix vector multiplication      |
!       |                y <--- OP*x                |
!       | The user should supply his/her own        |
!       | matrix vector multiplication routine here |
!       | that takes workd(ipntr(1)) as the input   |
!       | vector, and return the matrix vector      |
!       | product to workd(ipntr(2)).               | 
!       %-------------------------------------------%
        call cscmv(n, Hvalue, Hrow, Hpntrb, Hpntre, workd(ipntr(1):), workd(ipntr(2):))
    enddo
    
    if ( info < 0 ) then
!       %--------------------------%
!       | Error message, check the |
!       | documentation in ZNAUPD  |
!       %--------------------------%
        open( newunit = filerr, file='error_naupd.txt', action='write' )
        write(filerr, '(A,I0)') 'Error with _naupd, info = ', info
        write(filerr, '(A)')    'Check the documentation of _naupd'
        close(filerr)
        return
    endif
    
!   %-------------------------------------------%
!   | Post-Process using ZNEUPD .               |
!   |                                           |
!   | Computed eigenvalues may be extracted.    |
!   |                                           |
!   | Eigenvectors may also be computed now if  |
!   | desired.  (indicated by rvec = .true.)    |
!   %-------------------------------------------%
    call zneupd(rvec, 'A', select, d, v, ldv, sigma, &
                workev, bmat, n, which, nev, tol, resid, ncv, &
                v, ldv, iparam, ipntr, workd, workl, lworkl, &
                rwork, ierr)
!   %----------------------------------------------%
!   | Eigenvalues are returned in the one          |
!   | dimensional array D.  The corresponding      |
!   | eigenvectors are returned in the first NCONV |
!   | (=IPARAM(5)) columns of the two dimensional  | 
!   | array V if requested.  Otherwise, an         |
!   | orthogonal basis for the invariant subspace  |
!   | corresponding to the eigenvalues in D is     |
!   | returned in V.                               |
!   %----------------------------------------------%
    if ( ierr /= 0) then
!       %------------------------------------%
!       | Error condition:                   |
!       | Check the documentation of ZNEUPD. |
!       %------------------------------------%
        open( newunit = filerr, file='error_neupd.txt', action='write' )
        write(filerr, '(A,I0)') 'Error with _neupd, info = ', ierr
        write(filerr, '(A)')    'Check the documentation of _neupd.'
        close(filerr)
        return
    endif
    
    nconv = iparam(5)
    do j = 1, nconv
!       %---------------------------%
!       | Compute the residual norm |
!       |                           |
!       |   ||  A*x - lambda*x ||   |
!       |                           |
!       | for the NCONV accurately  |
!       | computed eigenvalues and  |
!       | eigenvectors.  (iparam(5) |
!       | indicates how many are    |
!       | accurate to the requested |
!       | tolerance)                |
!       %---------------------------%
        call cscmv(n, Hvalue, Hrow, Hpntrb, Hpntre, v(:,j), ax)
        call zaxpy (n, -d(j), v(1,j), 1, ax, 1)
        rd(j,1) = dble (d(j))
        rd(j,2) = aimag (d(j))
        rd(j,3) = dznrm2 (n, ax, 1)
        rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
    enddo
!   %-----------------------------%
!   | Display computed residuals. |
!   %-----------------------------%
    call dmout(filres, nconv, 3, rd, ncv, -6, 'Ritz values (Real, Imag) and relative residuals')
    
!   %-------------------------------------------%
!   | Print additional convergence information. |
!   %-------------------------------------------%
    if ( info == 1) then
        write(filconv, '(A)')
        write(filconv, '(A)') ' Maximum number of iterations reached.'
        write(filconv, '(A)')
    else if ( info == 3) then
        write(filconv, '(A)')
        write(filconv, '(A)') ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
        write(filconv, '(A)')
    end if
    
    write(filconv, '(A,I0)')      ' Size of the matrix :                           ', n
    write(filconv, '(A,I0)')      ' Number of Ritz values requested :              ', nev
    write(filconv, '(A,I0)')      ' Number of Arnoldi vectors generated (NCV) :    ', ncv
    write(filconv, '(2A)')        ' Portion of the spectrum :                      ', which
    write(filconv, '(A,I0)')      ' Number of converged Ritz values :              ', nconv
    write(filconv, '(A,I0)')      ' Number of Implicit Arnoldi update iterations : ', iparam(3)
    write(filconv, '(A,I0)')      ' Number of OP*x :                               ', iparam(9)
    write(filconv, '(A,1pG13.6)') ' Convergence criterion :                       ',  tol
    write(filconv, '(A)')
end subroutine eig

end module eigen
