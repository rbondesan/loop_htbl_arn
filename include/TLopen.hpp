// Class TL open. TL open boundary conditions (adjoint)

#ifndef TLOPEN_HPP
#define TLOPEN_HPP

#include <iostream>
#include <vector>
#include "hamiltonian.hpp"
#include "loop_htbl_states.hpp"
#include "parameters.hpp"
#include "common.hpp"
#include <cstdlib>
#include <complex>

using namespace std;

class TLopen : public Hamiltonian
{
public:
  TLopen (const Parameters& pars); 
  ~TLopen();
  
  int act(size_t cur_state, LoopHtblStates &states, 
	  vector< complex<double> > &tmp_nzval, vector<int> &tmp_irow, 
	  vector<int> &tmp_pcol);

};

#endif // TLOPEN_HPP
