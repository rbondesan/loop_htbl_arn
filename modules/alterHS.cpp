#include "alterHS.hpp"

//////////////////////////////////////////////////////////////
//                                                          //
//             Methods of class alterHS                     //
//                                                          //
//////////////////////////////////////////////////////////////

// ----------> Class:         alterHS
// ----------> Function name: alterHS
// ----------> Description:   Calls Hamiltonian constructor with pars
//
alterHS::alterHS(const Parameters& pars) : Hamiltonian (pars)
{ 
  cout << "Inside alterHS constructor" << endl; 
}

// ----------> Class:         alterHS
// ----------> Function name: ~alterHS
// ----------> Description:   Destructor
//
alterHS::~alterHS()
{
  cout << "Inside alterHS destructor" << endl;

}

// ----------> Class:         alterHS
// ----------> Function name: act
// ----------> Description:   act with the Hamiltonian, on states of position
//                            cur_state. Then update states and add the non
//                            zero element in tmp_nzval, tmp_irow and tmp_col
//
//
int alterHS::act(size_t cur_state, LoopHtblStates &states, 
		    vector< complex<double> > &tmp_nzval, vector<int> &tmp_irow, 
		    vector<int> &tmp_pcol)
{
  VecIntCoeff in_conf_coeff;
  VecIntCoeff out_conf_coeff;
  size_t start, end;
  size_t i,j;
  int first_pos_col;
  double cur_coupling;

  // Set tmp_col: it is the pointer to the values of tmp_nzval and
  // tmp_irow in column cur_state.
  // Hamiltonian::act acts on the given column cur_state
  tmp_pcol.push_back(tmp_nzval.size());
  first_pos_col = tmp_pcol.at( tmp_pcol.size()-1 );

  in_conf_coeff.conf = states.get_conf_at(cur_state);
  in_conf_coeff.coeff = 1;

  start = 0;
  end = 2*L; 
  for (i = start; i < end; i++)
    {
      for (j = i + 1; j < end; j++)
	{
	  //need type casting of size_t i,j
	  cur_coupling = 1./(
	    (2*L/M_PI * sin( (M_PI/(2*L))*((double) i- (double) j) ) )*
	    (2*L/M_PI * sin( (M_PI/(2*L))*((double) i- (double) j) ) ) );
	  // cout << "i = " << i << ", j = " << j << ", cur_coupling = " 
	  //       << cur_coupling << endl;
	  if (parity(i) != parity(j))
	    //i even, j odd or i even, j odd 
	    actTL(i, j, in_conf_coeff, out_conf_coeff,-cur_coupling);
	  else if (parity(i) == parity(j))
	    {   //i,j odd or i,j even
	      if (g != 0)
		actPerm(i, j, in_conf_coeff, out_conf_coeff,-cur_coupling*g);
	      else
		continue;
	    }
	  else
	    {
	      cerr << "in alterHS::act, wrong index parity, exiting" << endl;
	      exit(1);
	    }
#ifdef DEBUG      
	  for (size_t k = 0; k < n_conf; k++)    
	    cout << "in_conf[" << k << "] = " << (in_conf_coeff.conf).at(k)
		 << endl;
	  cout << "in coeff " << (in_conf_coeff).coeff << endl;
	  for (size_t k = 0; k < n_conf; k++)    
	    cout << "act op " << i << "," << j << " : out_conf[" << k << "] = " 
		 << (out_conf_coeff.conf).at(k) << endl;
	  cout << "coeff " << out_conf_coeff.coeff << endl;
#endif 
	  //Project onto the corresponding momentum eigenspace, 
	  //if not allowed range, just fill
	  project_momentum_and_fill(momentum, out_conf_coeff, 
				    states, tmp_nzval, tmp_irow, 
				    first_pos_col);
	}
    }

  // Sort tmp_irow and adjust corresponding nzval,
  // needed by the diagonalization routines
  quicksort_array_val(tmp_irow, tmp_nzval, tmp_pcol.at(cur_state), tmp_irow.size()-1);

  return 0;
}

// ----------> Class:         alterHS
// ----------> Function name: project_momentum_and_fill
// ----------> Description:   act on cur_conf_coeff with projector onto
//                            momentum sector, and fill into the (projected) 
//                            Hamiltonian and htbl the resulting states
//                            If momentum is not between 0 and L-1, just fill
//
int alterHS::project_momentum_and_fill(int momentum, VecIntCoeff cur_conf_coeff, 
				       LoopHtblStates &states, 
				       vector< complex<double> > &tmp_nzval, 
				       vector<int> &tmp_irow, 
				       int first_pos_col)
{
  size_t p;
  VecIntCoeff res_conf_coeff;
  complex<double> multiply;

  //skip if coeff = 0
  //If coeff not zero, then fill, otherwise do not fill, since
  //this state is simply not produced and acting on it later would be 
  //a mistake!
  if (cur_conf_coeff.coeff == 0.)
    return 0;
  //else:
  if (momentum > -1 && momentum < (int) L)
    {
      for (p = 0; p < L; p++)
	{
	  //act with translation of 2p sites coefficient omega^{2p}
	  //on cur_conf_coeff and produce res_conf_coeff
	  multiply = polar((double) 1/L, -2*M_PI/L*p*momentum);
	  act_translation(2*p, cur_conf_coeff, res_conf_coeff, multiply);
#ifdef DEBUG      
	  for (size_t j = 0; j < n_conf; j++)    
	    cout << "in_conf[" << j << "] = " << (cur_conf_coeff.conf).at(j)
		 << endl;
	  cout << "in coeff " << (cur_conf_coeff).coeff << endl;
	  for (size_t j = 0; j < n_conf; j++)    
	    cout << "act transl " << 2*p << ": out_conf[" << j << "] = " 
		 << (res_conf_coeff.conf).at(j) << endl;
	  cout << "coeff " << res_conf_coeff.coeff << endl;
#endif
	  //store res_conf_coeff with column given by first_pos_col - different
	  //tipically from that of cur_conf_coeff.
	  fill_nzval(states, (res_conf_coeff).conf, (res_conf_coeff).coeff, 
		     tmp_nzval, tmp_irow, first_pos_col);
	}
    }
  else //just fill if momentum is not in the allowed range.
    {
      fill_nzval(states, (cur_conf_coeff).conf, (cur_conf_coeff).coeff, 
		 tmp_nzval, tmp_irow, first_pos_col);
    }

  return 0;
}
