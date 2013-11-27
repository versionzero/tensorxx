#ifndef __COO3_H__
#define __COO3_H__

#include <algorithm>
#include <set>

/*
 * Compute B = A for COO tensor A, TMR(3) CSR tensor B
 *
 *
 * Input Arguments:
 *   I  n_slice    - number of (frontal) slices in A
 *   I  n_row      - number of rows in A
 *   I  n_col      - number of columns in A
 *   I  nnz        - number of nonzeros in A
 *   I  Ak[nnz(A)] - slice indices
 *   I  Ai[nnz(A)] - row indices
 *   I  Aj[nnz(A)] - column indices
 *   T  Ax[nnz(A)] - nonzeros
 * Output Arguments:
 *   I B_RO - row pointer
 *   I B_CO - column indices
 *   I B_KO - tube indices
 *   T B_V  - nonzeros
 *
 * Note:
 *   Output arrays B_RO, B_CO, B_KO, and B_V must be preallocated
 *
 * Note:
 *   Input: indices *are not* assumed to be ordered
 *
 *   Note: duplicate entries are carried over to the CSR represention
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col,n_slice))
 *
 */
template <class I, class T>
void coo_to_tmr_csr(const I n_slice,
		    const I n_row,
		    const I n_col,
		    const I nnz,
		    const I Ak[],
		    const I Ai[],
		    const I Aj[],
		    const T Ax[],
		    I       B_RO[],
		    I       B_CO[],
		    I       B_KO[],
		    T       B_V[])
{
  /* Compute number of non-zero entries per row of A */
  std::fill(B_RO, B_RO + n_row, 0);
  for (I n = 0; n < nnz; n++) {
    B_RO[Ai[n]]++;
  }
  
  /* Cumulative-sum the nnz per row to get B_RO[] */
  for (I i = 0, cumsum = 0; i < n_row; i++) {
    I temp  = B_RO[i];
    B_RO[i] = cumsum;
    cumsum += temp;
  }
  B_RO[n_row] = nnz; 
  
  /* Write Aj,Ak,Ax into B_CO,B_KO,B_V */
  for (I n = 0; n < nnz; n++) {
    I row      = Ai[n];
    I dest     = B_RO[row];
    B_CO[dest] = Aj[n];
    B_KO[dest] = Ak[n];
    B_V[dest]  = Ax[n];
    B_RO[row]++;
  }
  
  for (I i = 0, last = 0; i <= n_row; i++) {
    I temp  = B_RO[i];
    B_RO[i] = last;
    last    = temp;
  }
  
  /* Now B_RO,B_CO,B_V form a CSR representation (with possible
     duplicates) */
}

/*
 * Compute B = A for COO tensor A, ECSR tensor B
 *
 *
 * Input Arguments:
 *   I  n_row      - number of rows in A
 *   I  n_col      - number of columns in A
 *   I  n_slice     - number of tubes in A
 *   I  nnz        - number of nonzeros in A
 *   I  Ai[nnz(A)] - row indices
 *   I  Aj[nnz(A)] - column indices
 *   I  Ak[nnz(A)] - tube indices
 *   T  Ax[nnz(A)] - nonzeros
 * Output Arguments:
 *   I B_RO        - row pointer
 *   I B_CO        - column/tube indices
 *   T B_V         - nonzeros
 *
 * Note:
 *   Output arrays B_RO, B_CO and B_V must be preallocated
 *
 * Note:
 *   Input: indices *are not* assumed to be ordered
 *
 *   Note: duplicate entries are carried over to the ECSR represention
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col,n_slice))
 *
 */
template <class I, class T>
void coo_to_ecsr(const I n_row,
		 const I n_col,
		 const I n_slice,
		 const I nnz,
		 const I Ai[],
		 const I Aj[],
		 const I Ak[],
		 const T Ax[],
		 I       B_RO[],
		 I       B_CK[],
		 T       B_V[])
{
  /* Compute number of non-zero entries per row of A */
  std::fill(B_RO, B_RO + n_row, 0);
  for (I i = 0; i < nnz; i++) {
    B_RO[Ai[i]]++;
  }
  
  /* Cumulative-sum the nnz per row to get B_RO[] */
  for (I i = 0, cumsum = 0; i < n_row; i++) {
    I temp  = B_RO[i];
    B_RO[i] = cumsum;
    cumsum += temp;
  }
  B_RO[n_row] = nnz;
  
  /* Write Aj, Ak, Ax into B_CK, B_V */
  for (I i = 0; i < nnz; i++) {
    I row      = Ai[i];
    I dest     = B_RO[row];
    B_CK[dest] = (Aj[i] * n_slice) + Ak[i];
    B_V[dest]  = Ax[i];
    B_RO[row]++;
  }
  
  for (I i = 0, last = 0; i <= n_row; i++) {
    I temp  = B_RO[i];
    B_RO[i] = last;
    last    = temp;
  }
  
  /* Now B_RO, B_CK, B_V form a ECSR representation (with possible
     duplicates) */
}

/*
 * Compute B += A for COO matrix A, dense matrix B
 *
 * Input Arguments:
 *   I  n_row           - number of rows in A
 *   I  n_col           - number of columns in A
 *   I  nnz             - number of nonzeros in A
 *   I  Ai[nnz(A)]      - row indices
 *   I  Aj[nnz(A)]      - column indices
 *   T  Ax[nnz(A)]      - nonzeros 
 *   T  B_V[n_row*n_col] - dense matrix
 *
 */
template <class I, class T>
void coo_todense(const I n_row,
                 const I n_col,
                 const I nnz,
                 const I Ai[],
                 const I Aj[],
                 const T Ax[],
                       T B_V[],
		 int fortran)
{
  if (!fortran) {
    for(I i = 0; i < nnz; i++){
      B_V[ n_col * Ai[i] + Aj[i] ] += Ax[i];
    }
  } else {
    for(I i = 0; i < nnz; i++){
      B_V[ n_row * Aj[i] + Ai[i] ] += Ax[i];
    }
  }
}


/*
 * Compute Y += A*X for COO matrix A and dense vectors X,Y
 *
 *
 * Input Arguments:
 *   I  nnz           - number of nonzeros in A
 *   I  Ai[nnz]       - row indices
 *   I  Aj[nnz]       - column indices
 *   T  Ax[nnz]       - nonzero values
 *   T  Xx[n_col]     - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row]     - output vector
 *
 * Notes:
 *   Output array Yx must be preallocated
 *
 *   Complexity: Linear.  Specifically O(nnz(A))
 * 
 */
template <class I, class T>
void coo_matvec(const I nnz,
	            const I Ai[], 
	            const I Aj[], 
	            const T Ax[],
	            const T Xx[],
	                  T Yx[])
{
    for(I n = 0; n < nnz; n++){
        Yx[Ai[n]] += Ax[n] * Xx[Aj[n]];
    }
}

/*
 * Count the number of occupied diagonals in COO matrix A
 *
 * Input Arguments:
 *   I  nnz             - number of nonzeros in A
 *   I  Ai[nnz(A)]      - row indices
 *   I  Aj[nnz(A)]      - column indices
 *
 */
template <class I>
I coo_count_diagonals(const I nnz,
                      const I Ai[],
                      const I Aj[])
{
    std::set<I> diagonals;
    for(I n = 0; n < nnz; n++){
        diagonals.insert(Aj[n] - Ai[n]);
    }
    return diagonals.size();
}


#endif
