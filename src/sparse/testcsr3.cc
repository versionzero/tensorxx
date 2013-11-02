
#include <iostream>
#include "sparse.h"

int main(int argc, char **argv)
{
  int L, M, N, nnz;
  float *Ax, *Bx;
  int *Ai, *Aj, *Ak, *Bp, *Bj, *Bk;
  
  coo3_tocsr3(L, M, N, nnz,
	      Ai, Aj, Ak, Ax,
	      Bp, Bj, Bk, Bx);
  
  return 0;
}
