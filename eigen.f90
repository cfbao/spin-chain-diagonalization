module eigen
    implicit none
contains
! ==========================================================================
! This routine is an interface to (original) ARPACK. It computes selected 
! eigenvalues (and corresponding eigenvectors if requested) of a complex 
! matrix in sparse coordinate format.
! It is modified from the ARPACK example program ZNDRV1.F.
! ==========================================================================
! Routines called
!   znaupd          ARPACK reverse communication interface routine.
!   zneupd          ARPACK routine that returns Ritz values and (optionally)
!                   Ritz vectors (i.e. approx. eigenvalues and eigenvectors)
!   mkl_zcoogemv    Level 2 Sparse BLAS that computes matrix-vector product 
!                   of a sparse general matrix in the coordinate format
!   dlapy2          LAPACK routine to compute sqrt(x**2+y**2) carefully.
!   dznrm2          Level 1 BLAS that computes the norm of a complex vector.
!   zaxpy           Level 1 BLAS that computes y <- alpha*x+y.
! ==========================================================================
! Arguments
! --------------------------------------------------------------------------
! HVALUE        On input, store numeric value of nonzero matrix elements
! HROW          On input, store row index of nonzero matrix elements
! HCOLUMN       On input, store column index of nonzero matrix elements
! N             On input, specify dimension of the matrix
! NHELE         On input, specify number of nonzero matrix elements
! NEV           On input, specify number of eigenvalues to be computed
! NCV           On input, specify workspace size (i.e. 2nd dimension of V)
!                   Allowed range:            NEV+2 <= NCV <= N
!                   Recommended (in ARPACK):    NCV >= 2*NEV
!                   Tested:                NCV >= max{ 10, 3*NEV }
!                   Optimal size:            Problem dependant
!                   See APARCK documentation for full details.
! D             On output, store eigenvalues
! V             On input, array providing workspace, dimension (N,NCV).
!               On output, store eigenvectors in first NEV columns 
!               if requested.
! RVEC          On input,   .FALSE. to surpress computing eigenvectors
!                           .TRUE.  to request computing eigenvectors
! FILCONV       On input, specify file "pointer" to "convergence-k.txt"
! FILRES        On input, specify file "pointer" to "residuals-k.txt"
! --------------------------------------------------------------------------
! Local Parameters
! --------------------------------------------------------------------------
! WHICH     'LM' -> want the NEV eigenvalues of largest magnitude.
!           'SM' -> want the NEV eigenvalues of smallest magnitude.
!           'LR' -> want the NEV eigenvalues of largest real part.
!           'SR' -> want the NEV eigenvalues of smallest real part.
!           'LI' -> want the NEV eigenvalues of largest imaginary part.
!           'SI' -> want the NEV eigenvalues of smallest imaginary part.
! ==========================================================================
! Except from what is documented ABOVE this line, other parts of the routine 
! should not normally be changed. 
! If necessary, refer to full ARPACK documentation.
! ==========================================================================
! Local Large Arrays (ones that might concern memory use)
! --------------------------------------------------------------------------
! Complex(8):   AX(N)
!               WORKD(3*N)
!               WORKEV(3*N)
!               RESID(N)
! ==========================================================================

subroutine eig(Hvalue,Hrow,Hpntrb,Hpntre,nev, D,V, rvec, filconv, filres)
    use constants, only: zero, one
    implicit none
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
    character :: matdescra(6) = ['H', 'L', 'N', 'F', 'N', 'N']
    
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
        call mkl_zcscmv('N', n, n, one, matdescra, Hvalue, Hrow, Hpntrb, Hpntre, workd(ipntr(1)), zero, workd(ipntr(2)))
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
        call mkl_zcscmv('N', n, n, one, matdescra, Hvalue, Hrow, Hpntrb, Hpntre, v(1,j), zero, ax)
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