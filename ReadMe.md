# Fast exact diagonalization of translation-invariant spin chain

A Fortran program for fast exact diagonalization of translation-invariant spin chains.

Users only need to modify three `user_*.f90` files to use this program for their own purposes.

## TODO:
much more documentations coming soon

### Remarks:
* Relative residuals as written in residuals-k.txt should be ignored (unless eigenvectors are computed).
* As spin chain gets larger, array sizee may exceed $2^{31} - 1$, 
and ILP64 interface layer must be used.
Use makefile-i8 in that case.
