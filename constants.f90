! Copyright (C) 2012, 2014-2017 - Chenfeng Bao
!
! This program is free software; you can redistribute it and/or modify it 
! under the terms of the GNU General Public License; either version 3 of 
! the License, or (at your option) any later version.
! You should have received a copy of the GNU General Public License 
! along with this program; if not, see <http://www.gnu.org/licenses>.

module constants
    implicit none
    complex(8), parameter :: one = (1.0D0, 0.0D0), zero = (0.0D0, 0.0D0), onei = (0.0D0, 1.0D0)
    real(8), parameter :: pi = 3.1415926535897932D0
    character, parameter :: tabchar = achar(9)
end module constants
