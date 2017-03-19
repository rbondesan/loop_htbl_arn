#ifndef COMMON_H
#define COMMON_H

#define STR_VAL 0
#define STR_VAL_CONTRACT 100
#define STR_NUM 50

#define TLOPEN 0
#define TLPERIODIC 1
#define ALTERHS 2
#define ALTERCS 3 //alternating cirac-sierra

#include<vector>
#include<complex>

using namespace std;

typedef struct VecIntCoeff_
{
  vector<int> conf;
  complex<double> coeff;
} VecIntCoeff;

#endif
