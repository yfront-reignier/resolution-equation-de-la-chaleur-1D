/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

extern double dnrm2_(int *n, double *x, int *incx);

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the tridiagonal Poisson operator
  for(int i = 0;i < (*la) * (*lab); i++) {
    AB[i] = 0.0;
  }

  // TriL
  for(int i = 1; i < (*la); i++) {
    AB[(*kv) + i * (*lab)] = -1.0;
  }
  
  // Diag
  for(int i = 0; i < (*la); i++) {
    AB[(*kv + 1) + i * (*lab)] = 2.0;
  }

  // TriU
  for(int i = 0; i < (*la) - 1; i++) {
    AB[(*kv + 2) + i * (*lab)] = -1.0;
  }
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0
  for(int i = 0;i < (*la) * (*lab); i++) {
    AB[i] = 0.0;
  }

  for(int i = 0; i < (*la); i++) {
    AB[(*kv) + i * (*lab)] = 1.0;
  }
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  // TODO: Compute RHS vector
  for(int i = 0;i < (*la); i++) {
    RHS[i] = 0.0;
  }

  RHS[0] += *BC0;
  RHS[(*la) - 1] += *BC1;

}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D

  // T(x) = T(O) + x(T(1) - T(0))
  for(int i = 0; i < (*la); i ++) {
    EX_SOL[i] = (*BC0) + X[i] * (*BC1 - *BC0);
  }
  
}  

void set_grid_points_1D(double* x, int* la){
  // TODO: Generate uniformly spaced grid points in [0,1]
  double h = 1.0 / (*la + 1.0);

  for(int i = 0; i < *la; i++) {
    x[i] = (i+1) * h;
  }
}

double relative_forward_error(double* x, double* y, int* la){
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)
  int inc = 1;
  double norm_y, norm_diff;
  
  // ||y||
  norm_y = dnrm2_(la, y, &inc);
  if (norm_y < 1e-15) return 0.0;
  
  // ||x - y||
  double *diff_temp = (double*)malloc(*la * sizeof(double));
  for (int i = 0; i < *la; ++i) diff_temp[i] = x[i] - y[i];
  norm_diff = dnrm2_(la, diff_temp, &inc);
  free(diff_temp);
  
  return norm_diff / norm_y;
}

int indexABCol(int i, int j, int *lab){
  // TODO: Return the correct index formula for column-major band storage
  return i + j * (*lab);
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  dgbtrf_(la, n, kl, ku, AB, lab, ipiv, info);
  return *info;
}
