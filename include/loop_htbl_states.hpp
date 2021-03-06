//Class LoopHtblStates

#ifndef LOOP_HTBL_STATES_HPP
#define LOOP_HTBL_STATES_HPP

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>
#include "fun_def_states_htbl.hpp"
#include "parameters.hpp"
#include "common.hpp"

using namespace std;

class LoopHtblStates
{
public:
  LoopHtblStates();
  LoopHtblStates(const Parameters& pars); // init the states and insert them into the hash table
  ~LoopHtblStates();
  
  int init_states_and_add(const Parameters& pars);
  int check_and_add(vector<int> new_conf, int &pos_found);
  size_t get_n_states() const {return conf_vec_id.size() ;}; 
  IntVecId get_state_at(int cur_state) const {return conf_vec_id.at(cur_state) ;}; 
  vector<int> get_conf_at(int cur_state) const {return conf_vec_id.at(cur_state).conf ;}; 
  void print_table();
  int print_out_table(ofstream &outfile);

private:
  Myset htbl; //hash table type: unordered_set
  size_t n_conf;
  vector<IntVecId> conf_vec_id; // conf[cur_state] is the vector with
                                // values for coding the connectivity: 
                                // int a; a[0]=1;...;a[nconf]=2;
                                // conf.push_back(a);

};

#endif // LOOP_HTBL_STATES_HPP
