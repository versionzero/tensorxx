
#include <iostream>
#include "sparse.h"

void tensor_fread_coordinate(FILE *file, int *Ai, int *Aj, int *Ak, float *Ax)
{
  uint                 i, j, k;
  int                  l, m, n, q, nnz;
  int                  result;
  double               d;
  tensor_t             *tensor;
  tensor_storage_coordinate_t *storage;
  coordinate_tuple_t   *tuples;
  
  debug("tensor_fread_coordinate(0x%x)\n", file);
  
  if (0 != (result = mm_read_tensor_coordinate_size(file, &l, &m, &n, &nnz))) {
    die("Failed to read tensor dimensions (%d).\n", result);
  }
  
  debug("tensor_fread_coordinate: non-zero values: actual=%d.\n", nnz);
  debug("tensor_fread_coordinate: l=%d, m=%d, n=%d.\n", l, m, n);
  
  for (q = 0; q < nnz; ++q) {
    if (4 != (result = fscanf(file, "%u %u %u %lg\n", &k, &i, &j, &d))) {
      die("Failed to process line %d of the input stream (%d).\n", q, result);
    }
    tuples[q].index   = q;
    tuples[q].i       = i;
    tuples[q].j       = j;
    tuples[q].k       = k;
    tensor->values[q] = d;
  }
  
  return tensor;
}

int main(int argc, char **argv)
{
  int L, M, N, nnz;
  float *Ax, *B_V;
  int *Ai, *Aj, *Ak, *B_RO, *B_CO, *B_KO, *B_KO;
  
  coo_to_tmr_csr(M, N, L, nnz,
		 Ai, Aj, Ak, Ax,
		 Bp, Bj, Bk, Bx);
  
  coo_to_ecsr(M, N, L, nnz,
	      Ai, Aj, Ak, Ax,
	      Bp, Bj, Bk, Bx);
  
  return 0;
}
