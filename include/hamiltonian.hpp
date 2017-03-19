// Class Hamiltonian

#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

// TODO: pos in size_t !!!

#include <iostream>
#include <vector>
#include "loop_htbl_states.hpp"
#include "parameters.hpp"
#include "common.hpp"
#include <cstdlib>
#include <complex>

#define DIM_TABLE_ORIG 128

using namespace std;

class Hamiltonian
{
public:
  Hamiltonian(const Parameters& pars);
  ~Hamiltonian();

  // Define act as virtual and set it to zero. Implemented in every specific model
  virtual int act(size_t cur_state, LoopHtblStates &states, 
		  vector< complex<double> > &tmp_nzval, vector<int> &tmp_irow, 
		  vector<int> &tmp_pcol) = 0;

  protected:
  int fill_nzval(LoopHtblStates &states, vector<int> out_conf, complex<double> coeff,
		 vector< complex<double> > &tmp_nzval, vector<int> &tmp_irow, 
		 int first_pos_col);
  void quicksort_array_val(vector<int> &array, vector< complex<double> > &val, 
			   int startIndex, int endIndex);
  int split_array_val(vector<int> &array, vector< complex<double> > &val, int pivot, 
				 int startIndex, int endIndex);  
  template<typename T> void swap(T& left, T& right) {T tmp(left); left = right; right = tmp;};
  int actPerm(size_t index_i, size_t index_j, VecIntCoeff &in, VecIntCoeff &out, 
	      complex<double> multiply);
  int actTL(size_t index_i, size_t index_j, VecIntCoeff &in, VecIntCoeff &out,
	    complex<double> multiply);
  int act_translation(size_t num_sites, VecIntCoeff &in, VecIntCoeff &out, 
		      complex<double> multiply);
  int copy_el_comp_partner(vector<int> in_conf, vector<int> &out_conf, int value_i,
			   int value_j, size_t pos_i, size_t pos_j, size_t &p_pos_i, 
			   size_t &p_pos_j, int &small_value_available);
  int put_std_form(vector<int> &vec);
  void print_vec(vector<int> v){ 
    size_t ii; for (ii = 0; ii < v.size(); ii++)
		 cout<<"v["<<ii<<"]="<< v[ii] <<endl;};

  //Variables (protected)
  size_t L;
  size_t K;
  size_t n_conf;
  int momentum;
  double g;
  complex<double> fugacity;
  int model;

};

#endif // HAMILTONIAN_HPP
