#ifndef DIAG_EXACT_HPP
#define DIAG_EXACT_HPP

#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <complex>

using namespace std;

int diag_exact(int n, int nnz, complex<double> *nzval, int *irow, int *pcol);
int insert_in_mat(complex<double> *mat, int pos, complex<double> val, int size);
int print_result(int size, gsl_vector_complex *eval, gsl_matrix_complex *evec);
int diagonalize_non_symm(complex<double> *data, int size, gsl_vector_complex *eval,
                         gsl_matrix_complex *evec);


#endif
