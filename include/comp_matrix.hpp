#ifndef _COMP_MATRIX_HPP
#define _COMP_MATRIX_HPP

#include <vector>
#include <complex>
#include "parameters.hpp"
#include "hamiltonian.hpp"
#include "TLopen.hpp"
#include "TLperiodic.hpp"
#include "alterHS.hpp"
#include "alterCS.hpp"
#include <iostream>
#include <fstream>

using namespace std;

int comp_matrix(complex<double>* &nzval, int &nnz, int* &irow, int* &pcol, int &n,
		Parameters pars, bool print_out, ofstream &outfile);

int assign_values(vector< complex<double> > tmp_nzval, int &nnz, complex<double>* &nzval,  
		  vector<int> tmp_irow, int* &irow, int &n, int n_states,
		  vector<int> tmp_pcol, int* &pcol);

#endif
