//Class parameters

#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include "common.hpp"
#include <vector>
#include <complex>

using namespace std;

class Parameters
{
public:
  Parameters();
  Parameters(const Parameters& params);
  ~Parameters();
  
  //Read the parameters
  int read_parameters(string filename);
  int get_model() const {return model;};
  size_t get_L() const {return L;};
  size_t get_K() const {return K;};
  size_t get_n_conf() const {return 2*L; };
  size_t get_neval() const {return neval; };
  int get_momentum() const {return momentum;};
  double get_g() const {return g;};
  complex<double> get_fugacity() const {return fugacity; };

private:
  int model;
  size_t L;
  size_t K;
  int momentum; //momentum, if -1 means all in once 
  double g; //coupling g of the ALTERHS model (real double)
  complex<double> fugacity; //fugacity should be treated as a complex parameter,
  //same level of parameters in the Hamiltonian
  size_t neval; // Number of eigenvalues to be computed

};

#endif // PARAMETERS_HPP
