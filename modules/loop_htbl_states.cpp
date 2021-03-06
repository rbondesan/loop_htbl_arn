#include "loop_htbl_states.hpp"

//////////////////////////////////////////////////////////////
//                                                          //
//             Methods of class LoopHtblStates              //
//                                                          //
//////////////////////////////////////////////////////////////

// ----------> Class:         LoopHtblStates
// ----------> Function name: LoopHtblStates
// ----------> Description:   empty constructor
//                            
//
LoopHtblStates::LoopHtblStates()
{
  cout << "Inside LoopHtblStates empty constructor" << endl;

}

// ----------> Class:         LoopHtblStates
// ----------> Function name: LoopHtblStates
// ----------> Description:   Init htbl and states according
//                            to pars
//
LoopHtblStates::LoopHtblStates(const Parameters& pars)
{
  cout << "Inside LoopHtblStates constructor" << endl;

  n_conf = pars.get_n_conf();
  //Init conf_vec and n_states
  init_states_and_add(pars); 
}

// ----------> Class:         LoopHtblStates
// ----------> Function name: ~LoopHtblStates
// ----------> Description:   Destructor
//
LoopHtblStates::~LoopHtblStates()
{
  cout << "Inside LoopHtblStates destructor" << endl;

  // Dynamical memory is controlled by vector, no need to use delete
}

// ----------> Class:         LoopHtblStates
// ----------> Function name: init_states_and_add
// ----------> Description:   init the first state of model,
//                            L-K/2 link pattern in the bulk, K non contractible
//                            lines
//
//
int LoopHtblStates::init_states_and_add(const Parameters& pars)
{
  int model;
  vector<int> init_conf;
  IntVecId init_conf_id;
  size_t i;
  size_t K;
  int val;

  model = pars.get_model();
  K = pars.get_K();
  val = 0;

  if (pars.get_L()==0)
    {
      cerr << "L=0 not implemented. Exiting..." << endl;
      exit(1);
    }

  // TLopen
  if (model == TLOPEN || model == TLPERIODIC)
    {
      // Set the first K not contractible lines
      for (i = 0; i < K; i++)
	init_conf.push_back(STR_VAL);
      //	init_conf.push_back(-(i+1));	



      // Set the remaning n_conf - K to link patterns only juxtaposed
      //val--;
      for (i = K; i < n_conf; i++)  // n_conf is not even
	{ 
	  if ( ((i - K) % 2) == 0) // Increment every two.. 
	    val++;
	  init_conf.push_back(val);
	}
    }
  else if (model == ALTERHS || model == ALTERCS)
    {
      //Set the first K not contractible lines
      //treat them as distinguishable since they can be permuted
      //give them values from -1 to -K
      for (i = 0; i < K; i++)
	init_conf.push_back(-(i+1));
	//
	//init_conf.push_back(STR_VAL);

      // Set the remaning n_conf - K to link patterns only juxtaposed
      //val--;
      for (i = K; i < n_conf; i++)  // n_conf is not even
	{ 
	  if ( ((i - K) % 2) == 0) // Increment every two.. 
	    val++;
	  init_conf.push_back(val);
	}
    }
  else
    {
      cerr << "In LoopHtblStates::init_states_and_add , "
	   << "model not implemented. Exiting..." << endl;
      exit(1);
    }

  //Some trial values
  //   init_conf.push_back(53);
  //   init_conf.push_back(51);
  //   init_conf.push_back(52);
  //   init_conf.push_back(1);
  //   init_conf.push_back(1);
  //   init_conf.push_back(54);

  // Add init_conf to conf_vec and init n_states
  init_conf_id.conf = init_conf;
  init_conf_id.id = 0;
  conf_vec_id.push_back(init_conf_id);

  // and add this to the table
  Myset tmp_htbl(N_BUCK_FIXED,
		 my_hash_fun(),
		 my_match());
  htbl = tmp_htbl;
  htbl.insert(conf_vec_id.at( conf_vec_id.size() - 1 ));

  return 0;
}

// ----------> Class:         LoopHtblStates
// ----------> Function name: check_and_add
// ----------> Description:   check if new_conf is present using the hash table
//                            and if it is not the case, add the state
//
int LoopHtblStates::check_and_add(vector<int> new_conf, int &pos_found)
{
  IntVecId tmp_conf_id; 

  tmp_conf_id.conf = new_conf;
  tmp_conf_id.id = 0; // Fake value, match only on conf

  Myset::iterator it = htbl.find(tmp_conf_id);

  if ((it != htbl.end()) == true)
    {
      // Element is present
      pos_found = (*it).id;
    }
  else
    {
      // Element is not present, add it
      tmp_conf_id.id = conf_vec_id.size(); // Last but one plus one
      conf_vec_id.push_back( tmp_conf_id );
      // and add this to the table
      htbl.insert(conf_vec_id.at( conf_vec_id.size() - 1 ));

      // Update pos_found to the new position
      pos_found = tmp_conf_id.id;
    }

  return 0;
}

// ----------> Class:         LoopHtblStates
// ----------> Function name: print_table
// ----------> Description:   print the table of states
//
void LoopHtblStates::print_table()
{
  cout << "*** Table of states ****" << endl;
  for (Myset::const_iterator it = htbl.begin(); it != htbl.end(); ++it)
    {
      cout << "id:" << (*it).id << endl;
      for (size_t i = 0; i < n_conf; i++)
	cout << "conf[" << i << "] = " << (*it).conf.at(i) << endl;
    }
  cout << "***************************" << endl;
  
  return;
}
                           
// ----------> Class:         LoopHtblStates
// ----------> Function name: print_out_table
// ----------> Description:   print the table of states on outfile
//
int LoopHtblStates::print_out_table(ofstream &outfile)
{
  outfile << "#Table of states " << endl;
  for (Myset::const_iterator it = htbl.begin(); it != htbl.end(); ++it)
    {
      outfile << (*it).id << endl; //id
      for (size_t i = 0; i < n_conf; i++)
	{
	  if ((*it).conf.at(i) < 0)
	    outfile << 0 << endl; //if a string, put zero!
	  else
	    outfile << (*it).conf.at(i) << endl; //i-th configuration 
	}
    }
  
  return 0;
}
