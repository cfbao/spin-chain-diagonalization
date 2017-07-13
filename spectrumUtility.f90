module spectrumUtility
    implicit none
contains
subroutine sort_eigval(spectrumRaw, ncell, spectrumK, spectrumAll)
    real(8),intent(in) :: spectrumRaw(:, 0:)
    integer,intent(in) :: ncell
    real(8),allocatable,dimension(:,:),intent(out) :: spectrumK
    real(8),allocatable,dimension(:),  intent(out) :: spectrumAll
    
    integer n, k, nev, info
    nev = size(spectrumRaw, 1)
    allocate( spectrumAll(nev * ncell) )
    allocate( spectrumK(nev, -(ncell/2):ncell/2) )
    n = 1
    do k = 0, ncell/2
        spectrumAll(n:n+nev-1) = spectrumRaw(:, k)
        n = n + nev
        spectrumK(:, k) = spectrumRaw(:, k)
        if (k /= 0) then
            spectrumK(:, -k) = spectrumRaw(:, k)
        endif
        if ( k /= 0 .and. 2*k /= ncell ) then
            spectrumAll(n:n+nev-1) = spectrumRaw(:, k)
            n = n + nev
        endif
    enddo
    
    call dlasrt('I', size(spectrumAll), spectrumAll, info)
    do k = -(ncell/2), ncell/2
        call dlasrt('I', nev, spectrumK(:,k), info)
    enddo
end subroutine sort_eigval

subroutine write_eigval_k(spectrumK, outfile, ncell, eShift)
    use constants, only: tabchar
    integer,intent(in) :: ncell
    real(8),intent(in) :: spectrumK(:,-(ncell/2):)
    integer,intent(in) :: outfile
    real(8),intent(in) :: eShift
    
    integer k, i
    
    do k = -(ncell/2), ncell/2
        write(outfile, '(I4)', advance='no') k
        do i = 1, size(spectrumK, 1)
            write(outfile, '(A,1pG17.10)', advance='no') tabchar, spectrumK(i, k) + eShift
        enddo
        write(outfile,'(A)')
    enddo
    write(outfile,'(A)')
end subroutine write_eigval_k
end module spectrumUtility
