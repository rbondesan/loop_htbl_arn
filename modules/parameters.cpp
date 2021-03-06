#include "parameters.hpp"

//////////////////////////////////////////////////////////////
//                                                          //
//               Methods of class Parameters                //
//                                                          //
//////////////////////////////////////////////////////////////

// ----------> Class:         Parameters
// ----------> Function name: Parameters
// ----------> Description:   Empty constructor
//
Parameters::Parameters()
{
  cout << "Inside Parameters empty constructor" << endl;

  model = 0;
  L = 0;
  K = 0;  
  momentum = 0;
  g = 0;
  fugacity = 0;
  neval = 0;
}

// ----------> Class:         Parameters
// ----------> Function name: Parameters
// ----------> Description:   Copy constructor
//
Parameters::Parameters(const Parameters& params)
{
  cout << "Inside Parameters copy constructor" << endl;

  model = params.model;
  L = params.L;
  K = params.K;
  momentum = params.momentum;
  g = params.g;
  fugacity = params.fugacity;
  neval = params.neval;
}

// ----------> Class:         Parameters
// ----------> Function name: ~Parameters
// ----------> Description:   Destructor
//
Parameters::~Parameters()
{
  cout << "Inside Parameters destructor" << endl;
}

// ----------> Class:         Parameters
// ----------> Function name: read_data
// ----------> Description:   it opens the file with input and read
//
// TODO: Different parameters for different models!
//
int Parameters::read_parameters(string filename)
{
  ifstream infile;
  string buffer;

//  filename="data/" + filename;
  infile.open(filename.data());
  if(!infile)
    {
      cerr << "Inside Data::read_data. unable to open file " << filename
           << " for reading data. Exiting..." << endl;
      exit(1);
    }

  cout << "Inside Parameter::read_data. File: " << filename << 
    " opened for reading"  << "parameter" << endl;

  infile >> model;
  getline (infile, buffer);
  infile >> L;
  getline (infile, buffer);
  infile >> K;
  getline (infile, buffer);
  infile >> momentum;
  getline (infile, buffer);
  infile >> g;
  getline (infile, buffer);
  infile >> fugacity;
  getline (infile, buffer);
  infile >> neval;

  infile.close();

  return 0;
}
                           
