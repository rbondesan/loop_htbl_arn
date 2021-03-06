#include "TLopen.hpp"

//////////////////////////////////////////////////////////////
//                                                          //
//               Methods of class TLopen                    //
//                                                          //
//////////////////////////////////////////////////////////////

// ----------> Class:         TLopen
// ----------> Function name: TLopen
// ----------> Description:   Calls Hamiltonian constructor with pars
//
TLopen::TLopen(const Parameters& pars) : Hamiltonian (pars)
{ 
  cout << "Inside TLopen constructor" << endl; 
}

// ----------> Class:         TLopen
// ----------> Function name: ~TLopen
// ----------> Description:   Destructor
//
TLopen::~TLopen()
{
  cout << "Inside TLopen destructor" << endl;

}

// ----------> Class:         TLopen
// ----------> Function name: act
// ----------> Description:   act with the Hamiltonian, on states of position
//                            cur_state. Then update states and add the non
//                            zero element in tmp_nzval, tmp_irow and tmp_col
//
//
int TLopen::act(size_t cur_state, LoopHtblStates &states, 
		vector< complex<double> > &tmp_nzval, vector<int> &tmp_irow, 
		vector<int> &tmp_pcol)
{
  VecIntCoeff in_conf_coeff;
  VecIntCoeff out_conf_coeff;
  size_t start, end;
  size_t i;
  int first_pos_col;

  // Set tmp_col: it is the pointer to the values of tmp_nzval and
  // tmp_irow in column cur_state.
  // Hamiltonian::act acts on the given column cur_state
  tmp_pcol.push_back(tmp_nzval.size());
  first_pos_col = tmp_pcol.at( tmp_pcol.size()-1 );

  //If lbndles_symm symmetrize before acting with H
  in_conf_coeff.conf = states.get_conf_at(cur_state);
  in_conf_coeff.coeff = 1;

  start = 0;
  end = 2*L-1;
  for (i = start; i < end; i ++)
    {
      actTL(i, i+1, in_conf_coeff, out_conf_coeff,-1);
#ifdef DEBUG      
	  for (size_t j = 0; j < n_conf; j++)    
	    cout << "in_conf[" << j << "] = " << (in_conf_coeff.conf).at(j)
		 << endl;
	  cout << "in coeff " << (in_conf_coeff).coeff << endl;
	  for (size_t j = 0; j < n_conf; j++)    
	    cout << "act TL " << i << ": out_conf[" << j << "] = " 
		 << (out_conf_coeff.conf).at(j) << endl;
	  cout << "coeff " << out_conf_coeff.coeff << endl;
#endif
	  fill_nzval(states, (out_conf_coeff).conf, (out_conf_coeff).coeff, 
		     tmp_nzval, tmp_irow, first_pos_col);
    }

  // Sort tmp_irow and adjust corresponding nzval,
  // needed by the diagonalization routines
  quicksort_array_val(tmp_irow, tmp_nzval, tmp_pcol.at(cur_state), tmp_irow.size()-1);

  return 0;
}
