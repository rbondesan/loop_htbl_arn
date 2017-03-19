#ifndef DIAG_ARN_HPP
#define DIAG_ARN_HPP

#include "arlnsmat.h"
#include "arlsnsym.h"
#include "lnsymsol.h"

#include <string>
#include <iostream>
#include <fstream>

//using namespace std;

int diag_arn(int n, int nnz, double *nzval, int *irow, int *pcol, int neval,
	     bool print_out, std::ofstream &outfile);

#endif
