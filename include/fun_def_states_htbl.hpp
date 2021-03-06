#ifndef FUN_DEF_STATES_HTBL_HPP
#define FUN_DEF_STATES_HTBL_HPP

// Function for the hash table.
#include <iostream>
#include <vector>
//#include <tr1/unordered_set>
#include <unordered_set>
#include <cmath>
#include <cstdlib>

#define N_BUCK_FIXED 3539

using namespace std;

typedef struct IntVecId_
{
  vector<int> conf;
  int id;
} IntVecId;

// TODO: Improve the fact that we cannot set the number of buckets
// from the extern.

struct my_hash_fun
{
  size_t operator()(const IntVecId &el) const
  {
    long long int hash_key = 0;
    int K = 251; // A good choice of K
    size_t i;
    int ret;
    size_t n_conf;
    int n_buck;

    n_buck = N_BUCK_FIXED;   
    n_conf = el.conf.size();

    hash_key = abs(el.conf.at(n_conf - 1));
    for (i = 0; i < n_conf; i++)
      {
        hash_key = ( abs(el.conf.at(n_conf - 1 - i)) % n_buck
                     + (K % n_buck) * (hash_key % n_buck) % n_buck) % n_buck;
        if (hash_key < 0 || isnan(hash_key) || isinf(hash_key))
          {
            cerr << "Error! Problem with the computation of the hash key!"
              " hash_key " << hash_key << " n_bucketsisor %d " << n_buck << endl;
            exit(1);
          }
      }

    ret = (size_t)hash_key;

    return ret;
  }
};

// is equal function
struct my_match
{
  bool operator()(const IntVecId &el1, const IntVecId &el2) const
  {
    size_t i;
    bool iseq = true;

    for(i = 0; i < el1.conf.size(); i++)
      {
        if(el1.conf[i] != el2.conf[i])
          {
            iseq = false;
            break;
          }
      }

    return iseq;
  }
};

// The hash table type:
//typedef std::tr1::unordered_set<IntVecId, my_hash_fun, my_match> Myset;
typedef std::unordered_set<IntVecId, my_hash_fun, my_match> Myset;


#endif
