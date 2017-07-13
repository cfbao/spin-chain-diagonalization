! Copyright (C) 2012, 2014-2017 - Chenfeng Bao
!
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License; either version 3 of 
! the License, or (at your option) any later version.
! You should have received a copy of the GNU General Public License 
! along with this program; if not, see <http://www.gnu.org/licenses>.

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
