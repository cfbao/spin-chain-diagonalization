! ===================================================================================
! This program finds the low energy spectrum of a LBP-site-translation-invariant 
! spin-(D-1)/2 chain with NSITE = LBP*NBP sites. Periodic boundary condition is 
! assumed. The crystal momentum K takes value in {0, 1, ..., NBP-1}.
! Spectrum is assumed to be symmetric under k -> -k, hence only sectors with 
! K = 0, 1, ..., NBP/2 are actually computed.
! The Hamiltonian has one parameter 'G' that can be tuned at this level.
! -----------------------------------------------------------------------------------
! Default precision:  Double  {real(8) = 8 bytes,  complex(8) = 16 bytes}
! ===================================================================================
! Libraries used:
!    ARPACK-NG: Relevant files have already been included.
!    MKL:       Support functions (may be replaced, see Remarks 4)
!               LAPACK
!               BLAS Level 1, 2, 3
!               Sparse BLAS Level 2 
! -----------------------------------------------------------------------------------
! Routines used:
! -----------------------------------------------------------------------------------
! FIND_CYCLES: Find "building blocks" of momentum basis
! CONSTRUCTH:  Explicitly construct the Hamiltonian in each k-sector. 
!              Accept one parameter 'G'. Store H in sparse coordinate format.
!              Hamiltonian must be specified in this routine.
! EIG:         Solve the eigenvalue problem for one k-sector.
!              It is an interface to ARPACK-NG.
! DSECND:      MKL support function to get CPU time.
! ===================================================================================
! Input parameters
! -----------------------------------------------------------------------------------
! LBP:    Length of each translation invariant block
! NBP:    Number of translation invariant block.
! NSITE:  = LBP*NBP. Total number of sites of the system.
! D:      Dimension of local Hilbert space
! NTERM:  Number of basis vectors involved after applying the Hamiltonian on a basis 
!         vector. Roughly equals the number of nonzero matrix elements in a column.
!         usually of the same order as LBP*NBP
! G:      A parameter in the Hamiltonian
! NEV:    Number of eigenvalue (and eigenvectors) requested in each k-sector
! NCV:    Size of workspace array. Determine the size of EIGENVECTOR.
!         Allowed range:           NEV+2 <= NCV <= dim of matrix (NCONFIG, see below)
!         Recommended (in ARPACK): NCV >= 2*NEV
!         Tested:                  NCV >= max{ 20, 3*NEV }
!         See APARCK routines for full details.
! MEMMAX: Memory available on computer in unit of GB
! ===================================================================================
! Constant parameters (better left unchanged)
! -----------------------------------------------------------------------------------
! VALUE_VECTOR:   Set to .FALSE., suppress computation of eigenvectors.
!                 Current implementation does not support exporting eigenvectors.
! ===================================================================================
! Main variables
! -----------------------------------------------------------------------------------
! I,K,N:        Counters and occational usage
!               K reserved for crystal momentum
! NCONFIG:      Number of spin configurations (mod translations) 
!                     ~ 2**(2*NBP)/NBP
! EIGENVALUE:   Temporarily store eigenvalues computed in each k-sector.
! SPECTRUMRAW:  Store low energy spectrum as they're computed, minimum processing
! SPECTRUMALL:  Store the whole low energy spectrum in ascending order (GS first)
! SPECTRUMK:    Store the whole low energy spectrum in each k sector in ascending 
!               order (GS first)
! GS:           Ground state energy.
! NCONFIG_EST:  Estimated value of NCONFIG
! MEM:          Estimated memory usage in unit of GB
! FILSPEC       "Pointer" to spectrum.txt
! FILCONVK      "Pointer" to convergence-K.txt
! FILRESIK      "Pointer" to residuals-K.txt
! FILERR:       "Pointer" to error file(s)
! FILTIME:      "Pointer" to timer file
! -----------------------------------------------------------------------------------
! Large arrays (ones that may concern memory usage)
! -----------------------------------------------------------------------------------
! Arrays of size of NCONFIG (in module cycleUtility):
!    CONFIG:       store all (indices of) configurations mod translations 
!                  (i.e. a representative of each config.)
!    PERIOD:       Store the period of each spin configurations
!                  E.g. under 1-site translation,
!                  000000 has period 1, 001001 has period 3, 010000 has period 6
! Arrays of size of NHELE ~ NCONFIG*NTERM/2:
!    HVALUE:       Value of nonzero matrix elements
!                  (Usually the largest array in the entire program)
!    HROW:         Row index of nonzero matrix elements
!    HPNTRB:       Pointers to column beginnings
!    HPNTRE:       Pointers to column ends
! Arrays of size of NCONGIF*NCV:
!    EIGENVECTOR:  Workspace for solving eigenvalue problem.
!                  If requested, store eigenvectors on return.
!                  (Could be the largest array if NEV and NCV are large)
! ===================================================================================
! CAUTIONS and Remarks:
! -----------------------------------------------------------------------------------
! 1. The only integers with explicit INTEGER(8) data type are the indices of spin 
!    configurations in position basis (e.g. CONFIG(*)).
!    Related subroutines/functions (GET_ALPHA, GET_INDEX_X, GET_REP) must also use 
!    the correct data type.
! 2. HAMILTONIAN.F90 must be modified to implement a different Hamiltonian. 
!    Simplified instructions can be found in the documentation of HAMILTONIAN.F90.
!    Explanations of the technical details and mathematical background can be found 
!    in the note "Constructing Hamiltonian in k-Space".
! 3. Relative residuals as written in residuals-k.txt should be ignored unless
!    eigenvectors are computed.
! 4. This program can be modified to run with standalone LAPACK and BLAS without MKL,
!    as long as sparse matrix routines and timing function DSECND are replaced.
!    LAPACK 3.0 or before has its own DSECND; Fortran Intrinsic has CPU_TIME(SEC).
! 5. As system size gets larger, array sizee may exceed 2^31-1, and ILP64 interface
!    layer must be used (default INTERGER is 8-byte). Use makefile-i8 in that case.
! ===================================================================================
! Original program written in August 2012 at 
!                   Kavli Institute for Theoretical Physics.
! Generalized and documented in June 2014 at 
!                   Perimeter Institute for Theoretical Physics.
! F90 part modified to be Fortran 2008 standard conforming in September 2015 at
!                   Perimeter Institute for Theoretical Physics.
! F90 part modified to utilize more modern Fortran features in May 2016 at
!                   Perimeter Institute for Theoretical Physics.
! Chenfeng Bao
! ===================================================================================

program main
    use constants, only: tabchar
    use cycleUtility, only: find_cycles, config
    use Hamiltonian, only: constructH
    use eigen, only: eig
    use spectrumUtility, only: sort_eigval, write_eigval_k
    implicit none
    ! ================
    ! Input parameters
    ! ================
    integer,parameter :: nbp = 15, lbp = 1, nsite = lbp * nbp
    integer,parameter :: d = 2, nterm = nsite+1
    integer,parameter :: nev = 10, ncv = min(d**nsite, max(3*nev, 20))
    real(8),parameter :: g = -0.0_8, memmax = 1
    
    ! ===================
    ! Constant parameters
    ! ===================
    real(8),parameter :: degtol = 1.0D-6
    logical,parameter :: value_vector = .false.
    
    ! ===============
    ! Timer variables
    ! 1: Time before main loop
    ! 2: Time constructing Hamiltonian
    ! 3: Time solving eigen problem
    ! 4: Time outputting spectrum
    ! 5: Total time
    ! ===============
    real(8) timer(5), start(5), finish(5), dsecnd
    
    integer i, k, n, nconfig
    integer filspec, filconvK, filresiK, filerr, filtime
    real(8) nconfig_est, mem, gs
    
    ! ==================
    ! Allocatable arrays
    ! ==================
    integer,   allocatable,dimension(:)   :: Hrow, Hpntrb, Hpntre
    complex(8),allocatable,dimension(:)   :: Hvalue, eigenvalue
    complex(8),allocatable,dimension(:,:) :: eigenvector
    real(8),   allocatable,dimension(:)   :: spectrumAll
    real(8),   allocatable,dimension(:,:) :: spectrumRaw, spectrumK
    
    ! =================================================
    ! Pre-run checking: estimate array size & RAM usage
    ! If estimation exceeds capability, exit program.
    ! =================================================
    nconfig_est = (1.0_8*d)**nsite / nbp + (1.0_8*d)**(nsite/2) / nbp * 2 + 2
    mem = nconfig_est * ( 16*((nterm+1)/2.0_8 + ncv + 8) + 8*((nterm+1)/2.0_8 + 4) ) / 1000**3
    if( mem > memmax ) then
        open( newunit = filerr, file='not_enough_RAM.txt',action='write')
        write(filerr, '(A,1pG10.3,A)') 'Estimated memory usage :', mem,    ' GB'
        write(filerr, '(A,1pG10.3,A)') 'Available memory       :', memmax, ' GB'
        close(filerr)
        error stop 'Error: insufficient memory. See error file.'
    endif
    
    ! ============
    ! Timer starts
    ! ============
    timer = 0
    start(5) = dsecnd()
    
    ! ===================================================
    ! Preparation for constructing Hamiltonian.
    ! Relate spin position basis to momentum basis.
    ! Find explicitly all "cycles"#
    ! # "cycle" = spin configuration modulo translation
    ! ---------------------------------------------------
    ! 1. Count number of cycles, store in NCONFIG
    ! 2. Find indices/periods of representative of cycles
    !    store them in array CONFIG/PERIOD
    ! ===================================================
    start(1) = dsecnd()
    call find_cycles(nbp, lbp, d)
    finish(1) = dsecnd()
    timer(1) = timer(1) + ( finish(1) - start(1) )
    
    ! ===================================
    !           M A I N  L O O P
    ! Solve the system in each K-subspace
    ! -----------------------------------
    ! Throughout main loop, there are - 
    ! Three files open: 
    !    1. spectrum-K.txt
    !    2. convergence-K.txt
    !    3. residuals-K.txt
    ! One counter active: 
    !    N: # of total energy levels 
    ! ===================================
    open( newunit=filspec,  file='spectrum.txt',      action='write')
    open( newunit=filconvK, file='convergence-K.txt', action='write')
    open( newunit=filresiK, file='residuals-K.txt',   action='write')
    write(filspec, '(5A)')            'g', tabchar, 'nbp', tabchar, 'nsite'
    write(filspec, '(F7.4,A,I0,A,I0)') g,  tabchar,  nbp,  tabchar,  nsite
    write(filspec, '(A)')
    
    nconfig = size(config)
    n = nint(nconfig*(nterm+1)/2.0_8)
    allocate( Hvalue(n), Hrow(n), Hpntrb(nconfig), Hpntre(nconfig) )
    allocate( spectrumRaw(nev, 0:nbp/2) )
    main_loop: do k = 0, nBp/2
        ! =====================
        ! Construct Hamiltonian
        ! =====================
        start(2) = dsecnd()
        write(*, '(A,I0)') 'constructing Hamiltonian   k = ', k
        call constructH(Hvalue, Hrow, Hpntrb, Hpntre, nterm, k, g)
        finish(2) = dsecnd()
        timer(2) = timer(2) + ( finish(2) - start(2) )
        write(*, '(1pG0.3,A)') finish(2) - start(2), ' seconds'
        
        ! ===================
        ! Solve eigen problem
        ! ===================
        start(3) = dsecnd()
        write(*, '(A,I0)') 'solving eigenvalue problem k = ', k
        write(filconvK, '(A,I0)') 'k = ', k
        write(filresiK, '(A,I0)') 'k = ', k
        allocate( eigenvalue(ncv), eigenvector(nconfig, ncv) )
        call eig(Hvalue, Hrow, Hpntrb, Hpntre, nev, eigenvalue, eigenvector, &
                    value_vector, filconvK, filresiK)
        finish(3) = dsecnd()
        timer(3) = timer(3) + ( finish(3) - start(3) )
        write(*, '(1pG0.3,A)') finish(3) - start(3), ' seconds'
        
        ! ===============
        ! Export spectrum
        ! ===============
        start(4) = dsecnd()
        write(*, '(A,I0)') 'writing spectrum to file   k = ', k
        spectrumRaw(:, k) = eigenvalue(1:nev)
        write(filspec, '(I4)', advance='no') k
        do i = 1, nev
            write(filspec, '(A,1pG17.10)', advance='no') tabchar, spectrumRaw(i, k)
        enddo
        write(filspec, '(A)')
        finish(4) = dsecnd()
        timer(4) = timer(4) + ( finish(4) - start(4) )
        
        deallocate( eigenvalue, eigenvector )
    enddo main_loop
    write(filspec, '(A)')
    deallocate( Hvalue, Hrow, Hpntrb, Hpntre )
    close(filconvK)
    close(filresiK)
    ! ================
    ! End of MAIN LOOP
    ! ================
    ! Export spectrum
    ! ================
    start(4) = dsecnd()
    call sort_eigval(spectrumRaw, nbp, spectrumK, spectrumAll)
    gs = minval(spectrumAll)
    
    call write_eigval_k(spectrumK, filspec, nbp, -gs)
    call write_eigval_k(spectrumK, filspec, nbp, 0.0_8)
    do i = 1, size(spectrumAll)
        write(filspec, '(1pG17.10)') spectrumAll(i)
    enddo
    close(filspec)
    finish(4) = dsecnd()
    timer(4) = timer(4) + ( finish(4) - start(4) )
    
    ! =========================
    ! Timer ends
    ! Export timing information
    ! =========================
    finish(5) = dsecnd()
    timer(5) = timer(5) + ( finish(5) - start(5) )
    open( newunit = filtime, file='timer.txt', action='write' )
    write(filtime, '(A)') 'Time used:'
    write(filtime, '(A,1pG10.3,A)') 'Before main loop:         ', timer(1), ' seconds'
    write(filtime, '(A,1pG10.3,A)') 'Constructing Hamiltonian: ', timer(2), ' seconds'
    write(filtime, '(A,1pG10.3,A)') 'Solving eigen problem:    ', timer(3), ' seconds'
    write(filtime, '(A,1pG10.3,A)') 'Exporting spectrum:       ', timer(4), ' seconds'
    write(filtime, '(A,1pG10.3,A)') 'Total:                    ', timer(5), ' seconds'
    close(filtime)
end program main
