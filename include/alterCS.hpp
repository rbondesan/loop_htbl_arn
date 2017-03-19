// Class alternate Cirac-Sierra

#ifndef ALTERCS_HPP
#define ALTERCS_HPP

#include <iostream>
#include <vector>
#include "hamiltonian.hpp"
#include "loop_htbl_states.hpp"
#include "parameters.hpp"
#include "common.hpp"
#include <cstdlib>
#include <complex>

using namespace std;

class alterCS : public Hamiltonian
{
public:
  alterCS (const Parameters& pars); 
  ~alterCS();
  
  int act(size_t cur_state, LoopHtblStates &states, 
	  vector< complex<double> > &tmp_nzval, vector<int> &tmp_irow, 
	  vector<int> &tmp_pcol);
  int act2body(size_t index_i, size_t index_j, VecIntCoeff &in, VecIntCoeff &out,
	    complex<double> multiply);
  int act3body(size_t index_i, size_t index_j, size_t index_k,
	       VecIntCoeff &in, VecIntCoeff &out, complex<double> multiply);
  int project_momentum_and_fill(int momentum, VecIntCoeff cur_conf_coeff, 
				LoopHtblStates &states, 
				vector< complex<double> > &tmp_nzval, 
				vector<int> &tmp_irow, 
				int first_pos_col);
  size_t parity(size_t index){ return index % 2; };
  double mysign(size_t index){ if(index % 2 == 0) return 1.; else return -1.; };
  double mycot(double x){ return 1./tan(x); };

};

#endif // ALTERCS_HPP
