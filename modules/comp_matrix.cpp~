#include "comp_matrix.hpp"

/* A function to compute the non zero values of the matrix in CSC
   format.  non zero value are stored in the vector nzval of dimension
   nnz.  irow, pointer to vector that contains the row indices of non
   zero el stored in nzval.  irow, pointer to vector that contains the
   row indices of non zero el stored in nzval.
*/

int comp_matrix(complex<double>* &nzval, int &nnz, int* &irow, int* &pcol, int &n,
		Parameters pars, bool print_out, ofstream &outfile)
{
  // Use tmp values since the routines of arpack does not deal with
  // vectors type and we want to use this type since we don't know the
  // dimension of the vectors at the beginning
  vector< complex<double> > tmp_nzval; // Temporary values of nzval
  vector<int> tmp_irow;     // Temporary values of irow
  vector<int> tmp_pcol;     // Temporary values of pcol

  // State class. Constructor initializes hash table and states
  LoopHtblStates states(pars);
  size_t cur_state = 0;
  // Hamiltonian (abstract class)
  Hamiltonian * H;

  switch (pars.get_model()) 
    {
    case TLOPEN : 
      cout << "MODEL = TLOPEN" << endl;
      H = new TLopen(pars);
      break;
    case TLPERIODIC : 
      cout << "MODEL = TLperiodic" << endl;
      H = new TLperiodic(pars);
      break;
    case ALTERHS : 
      cout << "MODEL = alterHS" << endl;
      H = new alterHS(pars);
      break;
    case ALTERCS : 
      cout << "MODEL = alterCS" << endl;
      H = new alterCS(pars);
      break;
      //... 
    default : 
      cerr << "Error. No known model chosen in comp_matrix, exiting!"
	   << endl;
      exit(1);
    }

  // Cycle over states and fill tmp_nzval. Stopping criterium is that
  // n_states does not increment and we act on the last state
  while (cur_state != states.get_n_states()) 
    {
#ifdef DEBUG
      cout << "**** cur_state " << cur_state << " n_states "
	   << states.get_n_states() << " ****" << endl;
#endif
      H->act(cur_state, states, tmp_nzval, tmp_irow, tmp_pcol);
      cur_state++;
    }

#ifdef DEBUG
  states.print_table();
#endif

  //By convention, we define tmp_pcol[n+1] = nnz+1
  tmp_pcol.push_back(tmp_nzval.size());

  //delete the dynamically-allocated class
  //  delete[] H;

  // If here, we have compute the number of states in the model and tmp_nzval
  assign_values(tmp_nzval, nnz, nzval, tmp_irow, irow, n, states.get_n_states(),
		tmp_pcol, pcol);

  return 0;
}

// Assign the values to the return parameters
// The values in irow (and nzval) are sorted along the column
// Frees memory after each assignement
int assign_values(vector< complex<double> > tmp_nzval, int &nnz, complex<double>* &nzval,  
		  vector<int> tmp_irow, int* &irow, int &n, int n_states,
		  vector<int> tmp_pcol, int* &pcol)
{
  int i;

  nnz = tmp_nzval.size();
  n = n_states;
  nzval = new complex<double>[nnz];
  for (i = 0; i < nnz; i++)
    nzval[i] = tmp_nzval[i];
  tmp_nzval.clear();
  irow = new int[nnz];
  for (i = 0; i < nnz; i++)
    irow[i] = tmp_irow[i];
  tmp_irow.clear();
  pcol = new int[n+1];
  for (i = 0; i < n+1; i++)
    pcol[i] = tmp_pcol[i];
  tmp_pcol.clear();

  return 0;
}
