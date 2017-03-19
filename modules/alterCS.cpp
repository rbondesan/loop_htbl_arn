#include "alterCS.hpp"

//////////////////////////////////////////////////////////////
//                                                          //
//             Methods of class alterCS                     //
//                                                          //
//////////////////////////////////////////////////////////////

// ----------> Class:         alterCS
// ----------> Function name: alterCS
// ----------> Description:   Calls Hamiltonian constructor with pars
//
alterCS::alterCS(const Parameters& pars) : Hamiltonian (pars)
{ 
  cout << "Inside alterCS constructor" << endl; 
}

// ----------> Class:         alterCS
// ----------> Function name: ~alterCS
// ----------> Description:   Destructor
//
alterCS::~alterCS()
{
  cout << "Inside alterCS destructor" << endl;

}

// ----------> Class:         alterCS
// ----------> Function name: act
// ----------> Description:   act with the Hamiltonian, on states of position
//                            cur_state. Then update states and add the non
//                            zero element in tmp_nzval, tmp_irow and tmp_col
//
//
int alterCS::act(size_t cur_state, LoopHtblStates &states, 
		    vector< complex<double> > &tmp_nzval, vector<int> &tmp_irow, 
		    vector<int> &tmp_pcol)
{
  VecIntCoeff in_conf_coeff;
  VecIntCoeff out_conf_coeff;
  size_t start, end;
  size_t i,j,k;
  int first_pos_col;
  double gij,gijk;
  double di,dj,dk;
  double delta;
  double myconst;

  // Set tmp_col: it is the pointer to the values of tmp_nzval and
  // tmp_irow in column cur_state.
  // Hamiltonian::act acts on the given column cur_state
  tmp_pcol.push_back(tmp_nzval.size());
  first_pos_col = tmp_pcol.at( tmp_pcol.size()-1 );

  in_conf_coeff.conf = states.get_conf_at(cur_state);
  in_conf_coeff.coeff = 1;

  start = 0;
  end = 2*L; 
  delta = fugacity.real(); //no need to use complex
  for (i = start; i < end; i++)
    {
      for (j = i + 1; j < end; j++)
	{
	  di = (double) i;
	  dj = (double) j;
	  /////////////////
 	  // 2 body term //
 	  /////////////////
	  //need type casting of size_t i,j
	  //gij=0;
	  gij = -delta*(delta + 2) - 4*L*(delta + 1);
	  gij += -delta*delta*mysign(i)*mysign(j);
	  gij += (delta*delta*(1 + mysign(i)*mysign(j)) 
	   	     + 6*delta + 4)*1./(
					(sin( (M_PI/(2*L))*(di - dj) ) )*
					(sin( (M_PI/(2*L))*(di - dj) ) ) );
	  gij *= 1./2./(delta + 1);
	  // cout << "i = " << i << ", j = " << j << ", gij = " 
	  //       << gij << endl;
 	  act2body(i, j, in_conf_coeff, out_conf_coeff, gij);
	  //Project onto the corresponding momentum eigenspace, 
	  //if not allowed range, just fill
	  project_momentum_and_fill(momentum, out_conf_coeff, 
				    states, tmp_nzval, tmp_irow, 
				    first_pos_col);
 	  /////////////////
 	  // 3 body term //
 	  /////////////////
 	  for (k = j + 1; k < end; k++)
 	    {
	      dk = (double) k;
 	      gijk = -2*delta*(mysign(k) +
		  (mysign(k) - mysign(i))*
		  mycot((M_PI/(2*L))*(di - dj))*mycot((M_PI/(2*L))*(di - dk))
		  +
		  (mysign(k) - mysign(j))*
		  mycot((M_PI/(2*L))*(dj - di))*mycot((M_PI/(2*L))*(dj - dk))
		 );
	      gijk *= 1./4./(delta + 1);
	      // cout << "i = " << i << ", j = " << j 
	      //      << ", k " << k << ", gijk = " 
	      // 	   << gijk << endl;
	      //act always on the same in_conf_coeff.
	      act3body(i, j, k, in_conf_coeff, out_conf_coeff, gijk);
	      //Project onto the corresponding momentum eigenspace, 
	      //if not allowed range, just fill
	      project_momentum_and_fill(momentum, out_conf_coeff, 
					states, tmp_nzval, tmp_irow, 
					first_pos_col);
	      //Anitcommutator of operators:
	      act3body(j, i, k, in_conf_coeff, out_conf_coeff, gijk);
	      //Project onto the corresponding momentum eigenspace, 
	      //if not allowed range, just fill
	      project_momentum_and_fill(momentum, out_conf_coeff, 
					states, tmp_nzval, tmp_irow, 
					first_pos_col);
 	    }
	}
    }

  //add a constant which set gs energy to zero.
  myconst = (2*L-2)*(2*L)*(2*L)*(4*delta*delta + 7*delta + 2)/(24*(delta+1));
  myconst -= (2*L-2)*(2*L)*(delta*delta + delta - 1)/(6*(delta+1));
  out_conf_coeff.conf = in_conf_coeff.conf;
  out_conf_coeff.coeff = in_conf_coeff.coeff*myconst;
  project_momentum_and_fill(momentum, out_conf_coeff, 
  			    states, tmp_nzval, tmp_irow, 
  			    first_pos_col);

  // Sort tmp_irow and adjust corresponding nzval,
  // needed by the diagonalization routines
  quicksort_array_val(tmp_irow, tmp_nzval, tmp_pcol.at(cur_state), tmp_irow.size()-1);

  return 0;
}

// ----------> Class:         alterCS
// ----------> Function name: act2body
// ----------> Description:   act as a 2 body interaction. depending on the parity of the
//                            site, act with -TL or Perm. Note the minus sign, since
//                            2 body is the quad casismir of gl.
//                            multiply is the coupling constant
//
int alterCS::act2body(size_t index_i, size_t index_j, VecIntCoeff &in, VecIntCoeff &out,
		      complex<double> multiply)
{
  if (parity(index_i) != parity(index_j))
    {
      //i even, j odd or i even, j odd 
      //NOTE THE MINUS SIGN
      actTL(index_i, index_j, in, out, -multiply);
    }
  else if (parity(index_i) == parity(index_j))
    {   
      //i,j odd or i,j even
      actPerm(index_i, index_j, in, out, multiply);
    }
  else
    {
      cerr << "in alterCS::act2body, wrong index parity, exiting" << endl;
      exit(1);
    }

#ifdef DEBUG      
  for (size_t k = 0; k < n_conf; k++)    
    cout << "in_conf[" << k << "] = " << (in.conf).at(k)
	 << endl;
  cout << "in coeff " << (in).coeff << endl;
  for (size_t k = 0; k < n_conf; k++)    
    cout << "act op " << index_i << "," << index_j << " : out_conf[" << k << "] = " 
	 << (out.conf).at(k) << endl;
  cout << "coeff " << out.coeff << endl;
#endif

  return 0;
}

// ----------> Class:         alterCS
// ----------> Function name: act3body
// ----------> Description:   act as a 3 body interaction. depending on the parity of the
//                            site, many possibilities (8). Note the minus signs.
//                            multiply is the coupling constant
//
int alterCS::act3body(size_t index_i, size_t index_j, size_t index_k, 
		      VecIntCoeff &in, VecIntCoeff &out,
		      complex<double> multiply)
{
  //use to store the intermediate result
  //first of two diagrams. (diagram acting on a state produces only one other state.)
  VecIntCoeff tmp; 

  if (parity(index_i) == 0 && parity(index_j) == 0 && parity(index_k) == 0)
    {
      actPerm(index_j, index_k, in, tmp, multiply);
      actPerm(index_i, index_j, tmp, out, 1);
    }
  else if (parity(index_i) == 0 && parity(index_j) == 0 && parity(index_k) == 1)
    {
      actTL(index_j, index_k, in, tmp, -multiply);
      actPerm(index_i, index_j, tmp, out, 1);
    }
  else if (parity(index_i) == 0 && parity(index_j) == 1 && parity(index_k) == 0)
    {
      actTL(index_i, index_j, in, tmp, -multiply);
      actTL(index_j, index_k, tmp, out, 1);
    }
  else if (parity(index_i) == 1 && parity(index_j) == 0 && parity(index_k) == 0)
    {
      actPerm(index_j, index_k, in, tmp, -multiply);
      actTL(index_i, index_j, tmp, out, 1);
    }
  else if (parity(index_i) == 0 && parity(index_j) == 1 && parity(index_k) == 1)
    {
      actTL(index_i, index_j, in, tmp, multiply);
      actPerm(index_j, index_k, tmp, out, 1);
    }
  else if (parity(index_i) == 1 && parity(index_j) == 0 && parity(index_k) == 1)
    {
      actTL(index_j, index_k, in, tmp, multiply);
      actTL(index_i, index_j, tmp, out, 1);
    }
  else if (parity(index_i) == 1 && parity(index_j) == 1 && parity(index_k) == 0)
    {
      actPerm(index_i, index_j, in, tmp, multiply);
      actTL(index_j, index_k, tmp, out, 1);
    }
  else if (parity(index_i) == 1 && parity(index_j) == 1 && parity(index_k) == 1)
    {
      actPerm(index_i, index_j, in, tmp, -multiply);
      actPerm(index_j, index_k, tmp, out, 1);
    }
  else
    {
      cerr << "in alterCS::act3body, wrong index parity, exiting" << endl;
      exit(1);
    }

#ifdef DEBUG      
  for (size_t k = 0; k < n_conf; k++)    
    cout << "in_conf[" << k << "] = " << (in.conf).at(k)
	 << endl;
  cout << "in coeff " << in.coeff << endl;
  for (size_t k = 0; k < n_conf; k++)    
    cout << "act op " << index_i << "," << index_j << "," << index_k 
	 << " : out_conf[" << k << "] = " 
	 << (out.conf).at(k) << endl;
  cout << "coeff " << out.coeff << endl;
#endif

  return 0;
}

// ----------> Class:         alterCS
// ----------> Function name: project_momentum_and_fill
// ----------> Description:   act on cur_conf_coeff with projector onto
//                            momentum sector, and fill into the (projected) 
//                            Hamiltonian and htbl the resulting states
//                            If momentum is not between 0 and L-1, just fill
//
int alterCS::project_momentum_and_fill(int momentum, VecIntCoeff cur_conf_coeff, 
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
	  //Notice L and not 2L, since translation by two lattice spacings.
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


