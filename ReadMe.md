# Fast exact diagonalization of translation-invariant spin chain

### Remarks:
* Relative residuals as written in residuals-k.txt should be ignored (unless eigenvectors are computed).
* As system size gets larger, array sizee may exceed 2^31-1, 
and ILP64 interface layer must be used.
Use makefile-i8 in that case.
