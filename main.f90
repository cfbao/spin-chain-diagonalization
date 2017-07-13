! Copyright (C) 2012, 2014-2017 - Chenfeng Bao
!
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License; either version 3 of 
! the License, or (at your option) any later version.
! You should have received a copy of the GNU General Public License 
! along with this program; if not, see <http://www.gnu.org/licenses>.

! =================================================================================
! This program finds the low energy spectrum of a translation-invariant spin chain,
! assuming periodic boundary condition.
! ---------------------------------------------------------------------------------
! Libraries used:
!    ARPACK-NG: Relevant files already included
!    LAPACK
!    BLAS Level 1, 2, 3
!    MKL Sparse BLAS Level 2 (only used in user_cscmv.f90; may be substituted)
! =================================================================================
program main
    use user_parameters, only: ncell, lcell, nsite, d, nterm, nev, ncv, g, memmax
    use constants, only: tabchar
    use necklaces, only: find_orbits, count_necklaces, norbits
    use Hamiltonian, only: constructH
    use eigen, only: eig
    use spectrumUtility, only: sort_eigval, write_eigval_k
    implicit none
    ! ===================
    ! Constant parameters
    ! ===================
    real(8),parameter :: degtol = 1.0D-6
    logical,parameter :: value_vector = .false.     ! don't calculate eigenvectors
    
    ! ===============
    ! Timer variables
    ! 1: Time before main loop
    ! 2: Time constructing Hamiltonian
    ! 3: Time solving eigen problem
    ! 4: Time outputting spectrum
    ! 5: Total time
    ! ===============
    real(8) timer(5), start(5), finish(5)
    
    integer i, k, nelem
    integer filspec, filconvK, filresiK, filerr, filtime
    real(8) mem, gs
    
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
    norbits = count_necklaces(ncell, d**lcell)
    nelem = nint(norbits*(nterm+1)/2.0_8)               ! number of matrix elements for a Hermitian Hamiltonian
    mem = 16.0_8 * (nelem + norbits*ncv + ncv) &        ! COMPLEX: Hvalue, eigenvector, eigenvalue
        + 8.0_8 * nev * (ncell + ncell/2+1 + ncell+1) & ! REAL: spectrumAll, spectrumRaw, spectrumK
        + storage_size(i)/8.0_8 * &                     ! INTEGER:
         (nelem + 4*norbits + 2*d**nsite)               !   Hrow, Hpntrb, Hpntre, 4 arrays in module NECKLACES
    if( mem > memmax ) then
        open( newunit = filerr, file='not_enough_RAM.txt',action='write')
        write(filerr, '(A,1pG10.3,A)') 'Estimated memory usage :', mem,    ' Bytes'
        write(filerr, '(A,1pG10.3,A)') 'Available memory       :', memmax, ' Bytes'
        close(filerr)
        error stop 'Error: insufficient memory. See error file.'
    endif
    
    ! ============
    ! Timer starts
    ! ============
    timer = 0
    start(5) = perf_clock()
    
    ! ===================================================
    ! initialize module necklaces
    ! ===================================================
    start(1) = perf_clock()
    call find_orbits()
    finish(1) = perf_clock()
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
    ! ===================================
    open( newunit=filspec,  file='spectrum.txt',      action='write')
    open( newunit=filconvK, file='convergence-K.txt', action='write')
    open( newunit=filresiK, file='residuals-K.txt',   action='write')
    write(filspec, '(A)', advance='no') 'g = '
    write(filspec, *) g
    write(filspec, '("ncell =   ",I0)') ncell
    write(filspec, '("nsite = ",I0)') nsite
    write(filspec, '(A)')
    allocate( Hvalue(nelem), Hrow(nelem), Hpntrb(norbits), Hpntre(norbits) )
    allocate( spectrumRaw(nev, 0:ncell/2) )
    main_loop: do k = 0, ncell/2
        ! =====================
        ! Construct Hamiltonian
        ! =====================
        start(2) = perf_clock()
        write(*, '(A,I0)') 'constructing Hamiltonian   k = ', k
        call constructH(Hvalue, Hrow, Hpntrb, Hpntre, k)
        finish(2) = perf_clock()
        timer(2) = timer(2) + ( finish(2) - start(2) )
        write(*, '(1pG0.3,A)') finish(2) - start(2), ' seconds'
        
        ! ===================
        ! Solve eigen problem
        ! ===================
        start(3) = perf_clock()
        write(*, '(A,I0)') 'solving eigenvalue problem k = ', k
        write(filconvK, '(A,I0)') 'k = ', k
        write(filresiK, '(A,I0)') 'k = ', k
        allocate( eigenvalue(ncv), eigenvector(norbits, ncv) )
        call eig(Hvalue, Hrow, Hpntrb, Hpntre, nev, eigenvalue, eigenvector, &
                    value_vector, filconvK, filresiK)
        finish(3) = perf_clock()
        timer(3) = timer(3) + ( finish(3) - start(3) )
        write(*, '(1pG0.3,A)') finish(3) - start(3), ' seconds'
        
        ! ===============
        ! Export spectrum
        ! ===============
        start(4) = perf_clock()
        write(*, '(A,I0)') 'writing spectrum to file   k = ', k
        spectrumRaw(:, k) = eigenvalue(1:nev)
        write(filspec, '(I4)', advance='no') k
        do i = 1, nev
            write(filspec, '(A,1pG17.10)', advance='no') tabchar, spectrumRaw(i, k)
        enddo
        write(filspec, '(A)')
        finish(4) = perf_clock()
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
    start(4) = perf_clock()
    call sort_eigval(spectrumRaw, ncell, spectrumK, spectrumAll)
    gs = minval(spectrumAll)
    
    call write_eigval_k(spectrumK, filspec, ncell, -gs)
    call write_eigval_k(spectrumK, filspec, ncell, 0.0_8)
    do i = 1, size(spectrumAll)
        write(filspec, '(1pG17.10)') spectrumAll(i)
    enddo
    close(filspec)
    finish(4) = perf_clock()
    timer(4) = timer(4) + ( finish(4) - start(4) )
    
    ! =========================
    ! Timer ends
    ! Export timing information
    ! =========================
    finish(5) = perf_clock()
    timer(5) = timer(5) + ( finish(5) - start(5) )
    open( newunit = filtime, file='timer.txt', action='write' )
    write(filtime, '(A)') 'Time used:'
    write(filtime, '(A,1pG10.3,A)') 'Before main loop:         ', timer(1), ' seconds'
    write(filtime, '(A,1pG10.3,A)') 'Constructing Hamiltonian: ', timer(2), ' seconds'
    write(filtime, '(A,1pG10.3,A)') 'Solving eigen problem:    ', timer(3), ' seconds'
    write(filtime, '(A,1pG10.3,A)') 'Exporting spectrum:       ', timer(4), ' seconds'
    write(filtime, '(A,1pG10.3,A)') 'Total:                    ', timer(5), ' seconds'
    close(filtime)
contains

function perf_clock()
    ! timing function, a wrapper of system_clock
    integer(8) cnt, cnt_rate
    real(8) perf_clock
    intrinsic system_clock
    call system_clock(cnt, cnt_rate)
    perf_clock = 1.0_8 * cnt / cnt_rate
end function perf_clock

end program main
