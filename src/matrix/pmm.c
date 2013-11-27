
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#define THREAD_COUNT 2
#define SIZE         8

int **A, **B, **C, P[SIZE][SIZE][SIZE][SIZE];

int **matrix_malloc(int n)
{
  int i;
  int *vals, **temp;
  
  vals = (int*)  malloc(n*n*sizeof(int));
  temp = (int**) malloc(n*sizeof(int*));
  
  for (i = 0; i < n; ++i) {
    temp[i] = &(vals[i*n]);
  }
  
  return temp;
}

void access_pattern_init(int n)
{
  int i, j, k, m;
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      for (k = 0; k < n; ++k) {
	for (m = 0; m < n; ++m) {
	  P[i][j][k][m] = 0;
	}
      }
    }
  }
}

void access_pattern_print(int n)
{
  int i, j, k, m;
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      printf("(%d, %d) = [\n", i, j);
      for (k = 0; k < n; ++k) {
	printf("%d: ", k);
	for (m = 0; m < n; ++m) {
	  printf("%4d ", P[i][j][k][m]);
	}
	printf("\n");
      }
      printf("\n");
    }
    printf("]\n");
  }
}

void matrix_init(int first_row, int row_count, int n)
{
  int i, j;
  
  srand(n);
  
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; ++j) {
      if (i >= first_row && i < first_row + row_count) {
	A[i][j] = (int)(17.0 * rand()/(RAND_MAX + 1.0));
	B[i][j] = (int)(17.0 * rand()/(RAND_MAX + 1.0));
      }
    }
  }
}

void matrix_print(int **M, int n)
{
  int i,j;
  
  printf("0x%x = [\n", (unsigned int) M);
  for (i = 0; i < n; ++i) {
    printf("%d: ", i);
    for (j = 0; j < n; ++j) {
      printf("%4d ", M[i][j]);
    }
    printf("\n");
  }
  printf("]\n");
}

void matrix_multiply(int first_row, int row_count, int n)
{
  int i, j, k, v;
  
  for (i = first_row; i < first_row + row_count; ++i) {
    for (j = 0; j < n; ++j) {
      v = 0;
      for (k = 0; k < n; ++k) {
	v += A[i][k] * B[k][j];
	
	P[i][j][i][k]++;
	P[i][j][k][j]++;
      }
      C[i][j] = v;
    }
  }

}

void* work(void *id)
{
  int i, j;
  int n         = SIZE;
  int my_id     = *((int*) id);
  int row_count = n/THREAD_COUNT;
  int first_row = (int) my_id*row_count;
  
  matrix_init(first_row, row_count, n);
  access_pattern_init(n);
  
  matrix_multiply(first_row, row_count, n);
  
  return NULL;
}

int main()
{
  int	    i, j, k, m;
  int       n = SIZE;
  pthread_t threads[THREAD_COUNT];
  int	    thread_ids[THREAD_COUNT];
  
  A = matrix_malloc(n);
  B = matrix_malloc(n);
  C = matrix_malloc(n);
  
  for (i = 0; i < THREAD_COUNT; ++i) {
    thread_ids[i] = i;
    if (pthread_create(&(threads[i]), NULL, work,
		       (void*) &(thread_ids[i]))) {
      printf("Error creating thread %d\n", i);
      return 1;
    }
  }
  
  for (i = 0; i < THREAD_COUNT; i++) {
    pthread_join(threads[i], NULL);
  }
  
  matrix_print(A, n);
  matrix_print(B, n);
  matrix_print(C, n);
  access_pattern_print(n);
    
  return 0;
}
