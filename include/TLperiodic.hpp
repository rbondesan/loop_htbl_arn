// Class TL periodic. TL periodic boundary conditions

#ifndef TLPERIODIC_HPP
#define TLPERIODIC_HPP

#include <iostream>
#include <vector>
#include "hamiltonian.hpp"
#include "loop_htbl_states.hpp"
#include "parameters.hpp"
#include "common.hpp"
#include <cstdlib>
#include <complex>

using namespace std;

class TLperiodic : public Hamiltonian
{
public:
  TLperiodic (const Parameters& pars); 
  ~TLperiodic();
  
  int act(size_t cur_state, LoopHtblStates &states, 
		    vector< complex<double> > &tmp_nzval, vector<int> &tmp_irow, 
		    vector<int> &tmp_pcol);

  int project_momentum_and_fill(int momentum, VecIntCoeff cur_conf_coeff, 
				LoopHtblStates &states, 
				vector< complex<double> > &tmp_nzval, 
				vector<int> &tmp_irow, 
				int first_pos_col);

};

#endif // TLPERIODIC_HPP
