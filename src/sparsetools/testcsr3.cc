
#include <iostream>
#include "sparsetools.h"

int main(int argc, char **argv)
{
  int M, N, P, nnz;
  float *Ax, *Bx;
  int *Ai, *Aj, *Ak, *Bp, *Bj, *Bk;
  
  coo3_tocsr(N, M, P, nnz, 
	     Ai, Aj, Ak, Ax,
	     Bp, Bj, Bk, Bx);
  
  return 0;
}
