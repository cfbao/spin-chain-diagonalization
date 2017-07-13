module user_parameters
    implicit none
    integer,parameter :: ncell = 15             ! number of unit cells
    integer,parameter :: lcell = 1              ! number of sites per unit cell
    integer,parameter :: nsite = lcell * ncell  ! total number of sites
    integer,parameter :: d = 2                  ! Hilbert space dimension on one site
    integer,parameter :: nterm = nsite+1        ! number of terms in Hamiltonian (see docs)
    integer,parameter :: nev = 10               ! number of eigenvalues to be calculated in each k-sector
    integer,parameter :: ncv = min(d**nsite, max(3*nev, 20))    ! number of workspace vectors (see docs)
    real(8),parameter :: g = 0                  ! adjustable parameter(s) in Hamiltonian. may be scalar or array
    real(8),parameter :: memmax = 1024**3       ! available memory on machine (bytes)
end module user_parameters
