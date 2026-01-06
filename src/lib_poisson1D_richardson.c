/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

extern double dnrm2_(int *n, double *x, int *incx);
extern void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);
extern void dgbmv_(char *trans, int *m, int *n, int *kl, int *ku, 
                   double *alpha, double *a, int *lda, double *x, int *incx,
                   double *beta, double *y, int *incy);

void eig_poisson1D(double* eigval, int *la){
  // TODO: Compute all eigenvalues for the 1D Poisson operator
  int n = (*la);
  for(int i = 0; i < n; i++) {
    eigval[i] = 2 - 2 * cos((i * M_PI)/n+1);
  }
}

double eigmax_poisson1D(int *la){
  // TODO: Compute and return the maximum eigenvalue for the 1D Poisson operator
  return 2 - 2 * cos(((*la) * M_PI)/(*la) + 1);
}

double eigmin_poisson1D(int *la){
  // TODO: Compute and return the minimum eigenvalue for the 1D Poisson operator
  return 2 - 2 * cos((0 * M_PI)/(*la) + 1);
}

double richardson_alpha_opt(int *la){
  // TODO: Compute alpha_opt
  return 2 / (eigmin_poisson1D(la) + eigmax_poisson1D(la));
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat

  int inc = 1;
  double *temp_matrix = (double*)malloc(*la * sizeof(double));
  double *r = (double*)malloc(*la * sizeof(double));
  double norm_b, norm_r;
  double alpha_blas, beta_blas;

  norm_b = dnrm2_(la, RHS, &inc);

  for(*nbite = 0; *nbite < *maxit; (*nbite)++) {
    alpha_blas = 1.0;
    beta_blas = 0.0;
    dgbmv_("N", la, la, kl, ku, &alpha_blas, AB, lab, X, &inc, &beta_blas, temp_matrix, &inc);

    for(int i = 0; i < *la; i++) {
      r[i] = RHS[i] - temp_matrix[i];
    }

    norm_r = dnrm2_(la, r, &inc);
    resvec[*nbite] = norm_r/norm_b;

    if(resvec[*nbite] < *tol) {
      (*nbite)++;
      break;
    }

    daxpy_(la, alpha_rich, r, &inc, X, &inc);
  }
  free(temp_matrix);
  free(r);
}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A

  for(int i = 0; i < (*la)* (*lab); i++) {
    MB[i] = 0.0;
  }

  for(int i = 0; i < (*la); i++) {
    MB[(*kv + 1) + i * (*lab)] = AB[(*kv + 1) + i * (*lab)];
  }
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
    for(int i = 0; i < (*la)* (*lab); i++) {
    MB[i] = 0.0;
  }

  for(int i = 0; i < (*la); i++) {
    MB[(*kv + 1) + i * (*lab)] = AB[(*kv + 1) + i * (*lab)];
  }

  for(int i = 1; i < (*la); i++) {
    MB[(*kv + 2) + i * (*lab)] = AB[(*kv + 2) + i * (*lab)];
  }
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iterative method
  int inc = 1;
  int NRHS = 1;
  int info;
  double *temp_matrix = (double*)malloc(*la * sizeof(double));
  double *r = (double*)malloc(*la * sizeof(double));
  double *z = (double*)malloc(*la * sizeof(double));
  int *ipiv = (int*)malloc(*la * sizeof(int));
  double norm_b, norm_r;
  double alpha_blas, beta_blas;

  norm_b = dnrm2_(la, RHS, &inc);

  for(*nbite = 0; *nbite < *maxit; (*nbite)++) {
    alpha_blas = 1.0;
    beta_blas = 0.0;
    dgbmv_("N", la, la, kl, ku, &alpha_blas, AB, lab, X, &inc, &beta_blas, temp_matrix, &inc);

    for(int i = 0; i < *la; i++) {
      r[i] = RHS[i] - temp_matrix[i];
    }

    norm_r = dnrm2_(la, r, &inc);
    resvec[*nbite] = norm_r/norm_b;

    if(resvec[*nbite] < *tol) {
      (*nbite)++;
      break;
    }

    for(int i = 0; i < (*la); i++) {
      z[i] = r[i];
    }

    dgbtrf_(la, la, kl, ku, MB, lab, ipiv, &info);
    if(info == 0) {
      dgbtrs_("N", la, kl, ku, &NRHS, MB, lab, ipiv, z, la, &info);
    }

    daxpy_(la, &alpha_blas, z, &inc, X, &inc);
  }
  free(temp_matrix);
  free(r);
  free(z);
  free(ipiv);
}

