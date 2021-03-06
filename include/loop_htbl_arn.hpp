#ifndef LOOP_HTBL_ARN_HPP
#define LOOP_HTBL_ARN_HPP

//For Arnoldi complex:
#include "arcomp.h" //The "arcomplex" (complex) type definition, in arpack++/include
#include "arlnsmat.h" //The ARluNonSymMatrix class definition, in arpack++/include
#include "arlscomp.h" //The ARluCompStdEig class definition, in arpack++/include
#include "lcompsol.h"  //The Solution function, in arpack++/examples/matrices/complex/
//Old Arnoldi, for nonsymm problem
//#include "arlnsmat.h"
//#include "arlsnsym.h"
//#include "lnsymsol.h"

#include "comp_matrix.hpp"
#include "parameters.hpp"
#include "diag_exact.hpp"
//temporarely disabled
//#include "diag_arn.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <complex>

using namespace std;

#define ARNOLDI 0
#define EXACT 1

int get_from_cmdline(int argc, char **argv, string &filename, 
		     int &m_id, bool &print_out);

#endif
