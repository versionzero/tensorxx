#ifndef __COO3_H__
#define __COO3_H__

#include <algorithm>
#include <set>

/*
 * Compute B = A for COO tensor A, TMR(3) CSR tensor B
 *
 *
 * Input Arguments:
 *   I  n_row      - number of rows in A
 *   I  n_col      - number of columns in A
 *   I  n_tube     - number of tubes in A
 *   I  nnz        - number of nonzeros in A
 *   I  Ai[nnz(A)] - row indices
 *   I  Aj[nnz(A)] - column indices
 *   I  Ak[nnz(A)] - tube indices
 *   T  Ax[nnz(A)] - nonzeros
 * Output Arguments:
 *   I Bp  - row pointer
 *   I Bj  - column indices
 *   I Bk  - tube indices
 *   T Bx  - nonzeros
 *
 * Note:
 *   Output arrays Bp, Bj, Bk, and Bx must be preallocated
 *
 * Note: 
 *   Input: indices *are not* assumed to be ordered
 *           
 *   Note: duplicate entries are carried over to the CSR represention
 *
 *   Complexity: Linear.  Specifically O(nnz(A) + max(n_row,n_col,n_tube))
 * 
 */
template <class I, class T>
void coo3_to_csr3(const I n_tube,
		  const I n_row,
		  const I n_col,
		  const I nnz,
		  const I Ai[],
		  const I Aj[],
		  const I Ak[],
		  const T Ax[],
                        I Bp[],
                        I Bj[],
                        I Bk[],
                        T Bx[])
{
  /* Compute number of non-zero entries per row of A */
  std::fill(Bp, Bp + n_row, 0);
  for (I n = 0; n < nnz; n++) {
    Bp[Ai[n]]++;
  }
  
  /* Cumulative-sum the nnz per row to get Bp[] */
  for (I i = 0, cumsum = 0; i < n_row; i++) {
    I temp  = Bp[i];
    Bp[i]   = cumsum;
    cumsum += temp;
  }
  Bp[n_row] = nnz; 
  
  /* Write Aj,Ak,Ax into Bj,Bk,Bx */
  for (I n = 0; n < nnz; n++) {
    I row    = Ai[n];
    I dest   = Bp[row];
    Bj[dest] = Aj[n];
    Bk[dest] = Ak[n];
    Bx[dest] = Ax[n];
    Bp[row]++;
  }
  
  for (I i = 0, last = 0; i <= n_row; i++) {
    I temp = Bp[i];
    Bp[i]  = last;
    last   = temp;
  }
  
  /* Now Bp,Bj,Bx form a CSR representation (with possible
     duplicates) */
}

template<class I, class T>
void coo3_tocsc(const I n_row,
      	       const I n_col,
      	       const I nnz,
      	       const I Ai[],
      	       const I Aj[],
      	       const T Ax[],
      	             I Bp[],
      	             I Bi[],
      	             T Bx[])
{ coo3_tocsr<I,T>(n_col, n_row, nnz, Aj, Ai, Ax, Bp, Bi, Bx); }

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
 *   T  Bx[n_row*n_col] - dense matrix
 *
 */
template <class I, class T>
void coo3_todense(const I n_row,
                 const I n_col,
                 const I nnz,
                 const I Ai[],
                 const I Aj[],
                 const T Ax[],
                       T Bx[],
		 int fortran)
{
    if (!fortran) {
        for(I n = 0; n < nnz; n++){
            Bx[ n_col * Ai[n] + Aj[n] ] += Ax[n];
        }
    }
    else {
        for(I n = 0; n < nnz; n++){
            Bx[ n_row * Aj[n] + Ai[n] ] += Ax[n];
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
void coo3_matvec(const I nnz,
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
I coo3_count_diagonals(const I nnz,
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
