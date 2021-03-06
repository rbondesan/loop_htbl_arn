#include "TLperiodic.hpp"

//////////////////////////////////////////////////////////////
//                                                          //
//             Methods of class TLperiodic                  //
//                                                          //
//////////////////////////////////////////////////////////////

// ----------> Class:         TLperiodic
// ----------> Function name: TLperiodic
// ----------> Description:   Calls Hamiltonian constructor with pars
//
TLperiodic::TLperiodic(const Parameters& pars) : Hamiltonian (pars)
{ 
  cout << "Inside TLperiodic constructor" << endl; 
}

// ----------> Class:         TLperiodic
// ----------> Function name: ~TLperiodic
// ----------> Description:   Destructor
//
TLperiodic::~TLperiodic()
{
  cout << "Inside TLperiodic destructor" << endl;

}

// ----------> Class:         TLperiodic
// ----------> Function name: act
// ----------> Description:   act with the Hamiltonian, on states of position
//                            cur_state. Then update states and add the non
//                            zero element in tmp_nzval, tmp_irow and tmp_col
//
//
int TLperiodic::act(size_t cur_state, LoopHtblStates &states, 
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

  in_conf_coeff.conf = states.get_conf_at(cur_state);
  in_conf_coeff.coeff = 1;

  start = 0;
  end = 2*L; 
  for (i = start; i < end; i ++)
    {
      //Last element is TL of index 2*L-1 acting on 2L-1 and 0
      actTL(i, (i+1)%n_conf, in_conf_coeff, out_conf_coeff,-1);
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
      //If momentum, then project onto the corresponding eigenspace
      if (momentum > -1 && momentum < (int) L)
	{
	  //skip if coeff = 0... but nonetheless add it, as done in case no
	  //momentum
	  if (out_conf_coeff.coeff == 0.)
	    fill_nzval(states, (out_conf_coeff).conf, (out_conf_coeff).coeff, 
	  	       tmp_nzval, tmp_irow, first_pos_col);
	  else
	    project_momentum_and_fill(momentum, out_conf_coeff, 
				      states, 
				      tmp_nzval, 
				      tmp_irow, 
				      first_pos_col);
	}
      else //else, just fill as usual
	{
	  fill_nzval(states, (out_conf_coeff).conf, (out_conf_coeff).coeff, 
		     tmp_nzval, tmp_irow, first_pos_col);
	}
    }

  // Sort tmp_irow and adjust corresponding nzval,
  // needed by the diagonalization routines
  quicksort_array_val(tmp_irow, tmp_nzval, tmp_pcol.at(cur_state), tmp_irow.size()-1);

  return 0;
}

// ----------> Class:         TLperiodic
// ----------> Function name: project_momentum_and_fill
// ----------> Description:   act on cur_conf_coeff with projector onto
//                            momentum sector, and fill into the (projected) 
//                            Hamiltonian and htbl the resulting states
//
int TLperiodic::project_momentum_and_fill(int momentum, VecIntCoeff cur_conf_coeff, 
					  LoopHtblStates &states, 
					  vector< complex<double> > &tmp_nzval, 
					  vector<int> &tmp_irow, 
					  int first_pos_col)
{
  size_t p;
  VecIntCoeff res_conf_coeff;
  complex<double> multiply;

  if (momentum < 0 || momentum > (int) L-1)
    {
      cerr << "Wrong value of momentum, in TLperiodic::project_momentum_and_fill, " 
	   << "exiting. " << endl;

      exit(1);
    } //otherwise, continue  
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

  return 0;
}
