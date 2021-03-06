#include "hamiltonian.hpp"

//////////////////////////////////////////////////////////////
//                                                          //
//               Methods of class Hamiltonian               //
//                                                          //
//////////////////////////////////////////////////////////////

// ----------> Class:         Hamiltonian
// ----------> Function name: Hamiltonian
// ----------> Description:   
//
Hamiltonian::Hamiltonian(const Parameters& pars)
{
  cout << "Inside Hamiltonian constructor" << endl;

  L = pars.get_L();
  K = pars.get_K();
  n_conf = pars.get_n_conf();
  momentum = pars.get_momentum();
  g = pars.get_g();
  fugacity = pars.get_fugacity();
  model = pars.get_model();

}

// ----------> Class:         Hamiltonian
// ----------> Function name: ~Hamiltonian
// ----------> Description:   Destructor
//
Hamiltonian::~Hamiltonian()
{
  cout << "Inside Hamiltonian destructor" << endl;

}

// ----------> Class:         Hamiltonian
// ----------> Function name: fill_nzval
// ----------> Description:   check and eventually add out_conf to the table,
//                            and fill properly tmp_nzval and tmp_irow
//
// TODO : Improve this algorithm...
//
int Hamiltonian::fill_nzval(LoopHtblStates &states, vector<int> out_conf, 
			    complex<double> coeff, vector< complex<double> > &tmp_nzval, 
			    vector<int> &tmp_irow, int first_pos_col)
{
  int pos_out_conf;
  bool prev_found = false;
  size_t i;

  // Find the position of out conf in the table
  states.check_and_add(out_conf, pos_out_conf);
  
  // Now we know that H[pos_out_conf][cur_state] += coeff
  for (i = first_pos_col; i < tmp_irow.size(); i++)
    {
      if (tmp_irow[i] == pos_out_conf)
	{
	  tmp_nzval[i] += coeff;
	  prev_found = true;
	  break;
	}
      else
	prev_found = false;
    }
  
  if (prev_found == false)
    {
      if (coeff != 0. ) //Comparison between complex double numbers
	{
	  tmp_nzval.push_back(coeff);
	  tmp_irow.push_back(pos_out_conf);
	}
      else
	{
	  //Do nothing
	}
    }

  return 0;
}

// ----------> Class:         Hamiltonian
// ----------> Function name: quicksort_array_val
// ----------> Description:   sort ascending using quicksort array from startIndex
//                            endIndex, and adjust correspondigly val
//
void Hamiltonian::quicksort_array_val(vector<int> &array, vector< complex<double> > &val, 
				      int startIndex, int endIndex)
{
  //pivot element is the leftmost element
  int pivot;
  int splitPoint;

  //if they are equal, it means there is
  if(endIndex > startIndex)              
    //only one element and quicksort's job
    //here is finished
    { 
      pivot = array.at(startIndex);
      splitPoint = split_array_val(array, val, pivot, startIndex, endIndex);
      //SplitArray() returns the position where
      //pivot belongs to
      array.at(splitPoint) = pivot;
      quicksort_array_val(array, val, startIndex, splitPoint-1);   //Quick sort first half
      quicksort_array_val(array, val, splitPoint+1, endIndex);        //Quick sort second half
    }
}

// ----------> Class:         Hamiltonian
// ----------> Function name: split_array_val
// ----------> Description:   used by sort_irow_nzval 
//
int Hamiltonian::split_array_val(vector<int> &array, vector< complex<double> > &val, 
				 int pivot, int startIndex, int endIndex)
{
  int leftBoundary = startIndex;
  int rightBoundary = endIndex;

  //shuttle pivot until the boundaries meet
  while(leftBoundary < rightBoundary)             
    {
      //keep moving until a lesser element is found
      //or until the leftBoundary is reached
      while( pivot < array.at(rightBoundary) && rightBoundary > leftBoundary)      
        {
	  //move left
          rightBoundary--;                  
        }
      swap(array.at(leftBoundary), array.at(rightBoundary));
      swap(val.at(leftBoundary), val.at(rightBoundary));

      //keep moving until a greater or equal element is found
      //or until the rightBoundary is reached
      while( pivot >= array.at(leftBoundary) && leftBoundary < rightBoundary)          
        {
	  //move right
          leftBoundary++;                   
        }
      swap(array.at(rightBoundary), array.at(leftBoundary));
      swap(val.at(rightBoundary), val.at(leftBoundary));
    }
  //leftBoundary is the split point because
  return leftBoundary; 
  //the above while loop exits only when
  //leftBoundary and rightBoundary are equal      
}

// ----------> Class:         Hamiltonian
// ----------> Function name: actPerm
// ----------> Description:   Act with Permutation generator of index index on in_conf and 
//                            produce out_conf * coeff
//
int Hamiltonian::actPerm(size_t index_i, size_t index_j, VecIntCoeff &in, VecIntCoeff &out, 
			 complex<double> multiply)
{
  int value_ind_i;
  int value_ind_j;
  size_t k;
  
  out.coeff = multiply*in.coeff;
  // If out_conf contains value, erase them:
  (out.conf).erase((out.conf).begin(), (out.conf).end());
  
  value_ind_i = in.conf[index_i];
  value_ind_j = in.conf[index_j];

  for (k = 0; k < n_conf; k++)
    {
      if (k == index_i)
	(out.conf).push_back(value_ind_j);
      else if (k == index_j)
	(out.conf).push_back(value_ind_i);
      else
	(out.conf).push_back((in.conf).at(k));
    }

  put_std_form(out.conf);

  return 0;
}

// ----------> Class:         Hamiltonian
// ----------> Function name: actTL
// ----------> Description:   Act with TL generator of index index on in_conf and 
//                            produce out_conf * coeff, acts on sites i,j
//
int Hamiltonian::actTL(size_t index_i, size_t index_j, VecIntCoeff &in, 
		       VecIntCoeff &out, complex<double> multiply)
{
  //size_t index_j;
  //index_j=(index_i+1)%n_conf;
  size_t partner_index_i, partner_index_j;
  int in_conf_ind_i, in_conf_ind_j;
  int value;
  int small_value_available;

  out.coeff = multiply*in.coeff;
  // If out_conf contains value, erase them:
  (out.conf).erase((out.conf).begin(), (out.conf).end());

  in_conf_ind_i = in.conf[index_i];
  in_conf_ind_j = in.conf[index_j]; 

  if (in_conf_ind_i == in_conf_ind_j && (in_conf_ind_i != STR_VAL) 
      && (in_conf_ind_i != STR_VAL_CONTRACT))
    {
      //Here a link pattern |_| . Assigns 
      //same fugacity to contractible and not contractible loops.
      out.conf = in.conf;
      out.coeff *= fugacity;
    }
  else // in_conf_ind_i != in_conf_ind_j
    {
      copy_el_comp_partner(in.conf, out.conf, in_conf_ind_i, in_conf_ind_j, 
			   index_i, index_j, partner_index_i, partner_index_j, 
			   small_value_available);
      if (partner_index_i == index_i) //string
	{ 
	  if (partner_index_j == (index_j)) //Two strings
	    { 
	      if (in_conf_ind_i < 0 && in_conf_ind_j < 0)
		{
		  // Case acting on two distinguishable strings
		  out.coeff = 0;
		}
	      else if (in_conf_ind_i == STR_VAL_CONTRACT
		       && in_conf_ind_j  == STR_VAL_CONTRACT) 
		//Case in which two contractible strings contract them
		{
		  out.conf[index_j] = small_value_available;
		  out.conf[index_i] = small_value_available;		  		}
	      else if ((in_conf_ind_i == STR_VAL && in_conf_ind_j == STR_VAL_CONTRACT)
		       || (in_conf_ind_i == STR_VAL_CONTRACT && in_conf_ind_j == STR_VAL)) 
		//Case in which one string is contractible, the other not, ??
		{
		  cerr << "Not implemented this case, exit" << endl;
		  exit(1);
 
		}
	      else // Cases in which at least one of the two lines is a bulk string
		{ 
		  // Do nothing, as zero
		  out.coeff = 0;
		}
	    } // partner_index_j != (index_j), but partner_index_i == index_i
	  else // | |_|, |_|_| 
	    {
	      value = out.conf[index_j];
	      out.conf[partner_index_j] = out.conf[index_i];
	      out.conf[index_i] = value;
	      }
	}
      else // partner_index_i != index_i
	{
	  value = out.conf[index_i];
	  out.conf[partner_index_i] = out.conf[index_j];
	  out.conf[index_j] = value;
	}
      
      /* Modified, put in std form */
      put_std_form(out.conf);
    }
  
  return 0;
}

// ----------> Class:         Hamiltonian
// ----------> Function name: act_translation
// ----------> Description:   Act with translation of num_sites to the right
//                            on in_conf and produce out_conf * multiply
//
int Hamiltonian::act_translation(size_t num_sites, VecIntCoeff &in, VecIntCoeff &out, 
				 complex<double> multiply)
{
  size_t i;
  int dummy;
  size_t new_pos;

  out.coeff = multiply*in.coeff;
  // If out_conf contains value, erase them:
  (out.conf).erase((out.conf).begin(), (out.conf).end());
  
  for (i = 0; i < n_conf; i++)
    {
      dummy = (n_conf - num_sites);
      new_pos = (i + (size_t) dummy) % n_conf;
      // cout << "2*p = " << num_sites << ", i = " << i <<  
      // 	", (i - 2*p) % n_conf " << new_pos << endl;
      (out.conf).push_back((in.conf).at(new_pos));
    }

  put_std_form(out.conf);

  return 0;
}

// OLD version, acts on index, index+1
// // ----------> Class:         Hamiltonian
// // ----------> Function name: actTL
// // ----------> Description:   Act with TL generator of index index on in_conf and 
// //                            produce out_conf * coeff
// //
// int Hamiltonian::actTL(size_t index, VecIntCoeff &in, VecIntCoeff &out, double multiply)
// {
//   size_t index_p1;
//   size_t partner_index, partner_index_p1;
//   int in_conf_ind, in_conf_ind_p1;
//   int value;
//   int small_value_available;

//   out.coeff = multiply*in.coeff;
//   //If index=n_conf-1, then index+1 has to be 0
//   index_p1 = (index + 1 )% n_conf;
//   // If out_conf contains value, erase them:
//   (out.conf).erase((out.conf).begin(), (out.conf).end());

//   in_conf_ind = in.conf[index];
//   in_conf_ind_p1 = in.conf[index_p1]; 

//   if (in_conf_ind == in_conf_ind_p1 && (in_conf_ind != STR_VAL) 
//       && (in_conf_ind != STR_VAL_CONTRACT))
//     {
//       //Here a link pattern |_| . Assigns 
//       //same fugacity to contractible and not contractible loops.
//       out.conf = in.conf;
//       out.coeff *= fugacity;
//     }
//   else // in_conf_ind != in_conf_ind_p1 
//     {
//       copy_el_comp_partner(in.conf, out.conf, in_conf_ind, in_conf_ind_p1, 
// 			   index, partner_index, partner_index_p1, small_value_available);
//       if (partner_index == index) //string
// 	{ 
// 	  if (partner_index_p1 == (index_p1)) //Two strings
// 	    { 
// 	      if (in_conf_ind < 0 && in_conf_ind_p1 < 0)
// 		{
// 		  // Case acting on two boundary strings
// 		}
// 	      else if (in_conf_ind == STR_VAL_CONTRACT
// 		       && in_conf_ind_p1  == STR_VAL_CONTRACT) 
// 		//Case in which two contractible strings contract them
// 		{
// 		  out.conf[index_p1] = small_value_available;
// 		  out.conf[index] = small_value_available;		  		}
// 	      else if ((in_conf_ind == STR_VAL && in_conf_ind_p1 == STR_VAL_CONTRACT)
// 		       || (in_conf_ind == STR_VAL_CONTRACT && in_conf_ind_p1 == STR_VAL)) 
// 		//Case in which one string is contractible, the other not, ??
// 		{
// 		  cerr << "Not implemented this case, exit" << endl;
// 		  exit(1);
 
// 		}
// 	      else // Cases in which at least one of the two lines is a bulk string
// 		{ 
// 		  // Do nothing, as zero
// 		  out.coeff = 0;
// 		}
// 	    } // partner_index_p1 != (index_p1), but partner_index == index 
// 	  else // | |_|, |_|_| 
// 	    {
// 	      value = out.conf[index_p1];
// 	      out.conf[partner_index_p1] = out.conf[index];
// 	      out.conf[index] = value;
// 	      }
// 	}
//       else // partner_index != index 
// 	{
// 	  value = out.conf[index];
// 	  out.conf[partner_index] = out.conf[index_p1];
// 	  out.conf[index_p1] = value;
// 	}
      
//       /* Modified, put in std form */
//       put_std_form(out.conf);
//     }
  
//   return 0;
// }

// ----------> Class:         Hamiltonian
// ----------> Function name: copy_el_comp_partner
// ----------> Description:   copy and computes partner (element that form an arc with the 
//                            element that has value = value in position pos. If a string,
//                            return the same position ): LINEAR TIME
//                            Know already that it is not in pos_p1, find also that partner.
//
int Hamiltonian::copy_el_comp_partner(vector<int> in_conf, vector<int> &out_conf, int value_i,
				      int value_j, size_t pos_i, size_t pos_j, size_t &p_pos_i, 
				      size_t &p_pos_j, int &small_value_available)
{
  size_t k;
  int found_pos_i = 0;
  int found_pos_j = 0;
  int cur_value;
  
  // First set p_pos_i and p_pos_j equal to pos_i and pos_j, so that if we do not find any
  // partner, this means that are strings 
  p_pos_i = pos_i;
  p_pos_j = pos_j;
  
  // Small value available 
  small_value_available = 1;

  // If a string in the bulk give to the partner, the same position 
  // Need this if we code the strings with the same STR_VAL value!!
  if (value_i == STR_VAL || value_i == STR_VAL_CONTRACT)
    found_pos_i = 1;
  if (value_j == STR_VAL || value_j == STR_VAL_CONTRACT)
    found_pos_j = 1;
  
  // Find the partners, copy el 
  for (k = 0; k < in_conf.size(); k++)
    { 
      cur_value = in_conf[k];
      out_conf.push_back(cur_value);      

      // TODO: ????
      if (cur_value >= small_value_available && cur_value < STR_NUM)
	small_value_available = cur_value + 1;

      if (found_pos_i == 0)
	{ 
	  if (cur_value == value_i && k != pos_i)
	    { 
	      p_pos_i = k;
	      found_pos_i = 1; //Found partner of element in pos_i
	    }
	}
      if (found_pos_j == 0)
	{
	  if (cur_value == value_j && k != pos_j)
	    {
	      p_pos_j = k;
	      found_pos_j = 1; //Found partner of element in pos_j
	    }
	}
    }
  
  return 0;
}

// ----------> Class:         Hamiltonian
// ----------> Function name: put_std_form
// ----------> Description:   Put vec in std form, increasing numbers codifying
//                            the link pattern, while do not touch the noncontractible
//                            lines, since their value ordered means their position.
//
int Hamiltonian::put_std_form(vector<int> &vec)
{
  int i;
  int dim;
  vector<int> table;
  int num_to_assign;
  
  dim = vec.size();

  if (dim > DIM_TABLE_ORIG/2) // ???? 
    {
      cerr << "Error in Hamiltonian::put_std_form! Exiting" << endl;
      exit(1);
    }
   
  /* Initialize table to zero */
  for (i = 0; i < DIM_TABLE_ORIG; i++)
    table.push_back(0);

  num_to_assign = 0;
  for (i = 0; i < dim; i++)
    {
      if (vec[i] == STR_VAL || vec[i] == STR_VAL_CONTRACT)
	{
	  /* Do nothing, strings in the bulk */
	}
      else if (vec[i] < 0)
	{
	  /* Strings on the boundary, do nothing */
	}
      else if (vec[i] > STR_NUM)
	{
	  /* Strings in the bulk distinguishable */
	}
      else
	{
	  if (table[vec[i]] != 0)
	    vec[i] = table[vec[i]];
	  else
	    {
	      num_to_assign++;
	      table[vec[i]] = num_to_assign;
	      vec[i] = num_to_assign;
	    }
	}
    }

  return 0;
}

