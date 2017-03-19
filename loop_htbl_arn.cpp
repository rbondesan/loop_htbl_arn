// loop_htbl_arn:
// diagonalize the hamiltonian in the loop representation
// using the Arnoldi method, accessed through the
// libraries arpack++, alternatively use GLS

#include "loop_htbl_arn.hpp"

int main(int argc, char **argv)
{
  string filename;
  Parameters cur_params; // Parameters of the problem
  int     n;             // Dimension of the problem.
  int     nnz;           // Number of nonzero elements in the matrix A.
  int*    irow;          // pointer to an array that stores the row
                         // indices of the nonzeros in A.
  int*    pcol;          // pointer to an array of pointers to the
                         // beginning of each column of A in vector A.
  complex<double>* nzval;// pointer to an array that stores the
                         // nonzero elements of A. (complex double)
  int method_id;         // Index of the method to be used
  //Not used now:
  bool print_out=false;  // Print a lot of informations on states on output file
  ofstream outfile;      //File to print to if print_out is true
  string outname;        //Name of file to print to if print_out is true
  stringstream ss (stringstream::in | stringstream::out); //String to save file name
  //End not used

  //Get input from command line
  get_from_cmdline(argc, argv, filename, method_id, print_out);
  //initialize cur_params, element of class parameter
  cur_params.read_parameters(filename);
  // get the non zero element of the matrix in nzval.
  // Allocate nzval, irow, pcol
  cout << " Starting the computation of the matrix " << endl;
  cout << " L " << cur_params.get_L() << " K " << cur_params.get_K()
       << " momentum " << cur_params.get_momentum() 
       << " fugacity " << cur_params.get_fugacity() 
       << " neval " << cur_params.get_neval() << endl;
  comp_matrix(nzval, nnz, irow, pcol, n, cur_params, print_out, outfile);
  cout << " Computation of the matrix ended: " << endl;
  cout << " n = " << n << " nnz " << nnz << endl;
  cout << endl;
  cout << " ******************************************* " << endl;
  cout << endl;

  //Start diagonalization
  cout << " Starting diagonalization of the matrix " << endl;
  if (method_id == ARNOLDI)
    {
      //cout << " Arnoldi Method temporarely disabled" << endl;
      cout << " Arnoldi Method for complex matrices " << endl;
      //////////////////////////
      ARluNonSymMatrix<arcomplex<double>, double> matrix(n, nnz, nzval, irow, pcol);
      //      ARluNonSymMatrix<double, double> matrix(n, nnz, nzval, irow, pcol);
      // Defining what we need: neval eigenvectors matrix A with largest
      // magnitude (LM)
      //      ARluNonSymStdEig<double> dprob(cur_params.get_neval(), matrix);
      ARluCompStdEig<double> dprob(cur_params.get_neval(), matrix);
      dprob.ChangeWhich("SR"); 
      //dprob.ChangeWhich("LR"); // LARGEST REAL
      // Turns on the debug mode
      dprob.Trace();
      // Printing solution.
      dprob.FindEigenvalues();
      Solution(matrix, dprob);      
      /////////////////////////
      //OLD
      //diag_arn(n, nnz, nzval, irow, pcol, cur_params.get_neval(),
      //        print_out, outfile);
    }
  else if (method_id == EXACT)
    {
      //      cout << " Exact Diagonalization temporarely disabled" << endl;
      //Diag exact does not support print out option...
      diag_exact(n, nnz, nzval, irow, pcol);
    }
  else 
    {
      cerr << "Please specify a method by -mMethod, Method=0,1! Exiting... " << endl;
      exit(1);
    }

  cout << " Diagonalization of the matrix ended " << endl;

  if (print_out == true)
    {
      outfile.close();
      cout << "Wrote on out file " << outname << endl;
    }

  // Delete the array previouly allocated
  delete [] nzval;
  delete [] irow;
  delete [] pcol;
  
  return 0;
}

// TODO: Use getopt...
int get_from_cmdline(int argc, char **argv, string &filename, int &m_id,
		     bool &print_out)
{
  stringstream ss(stringstream::in | stringstream::out);

  // Get the options
  if (argc == 1)
    {
      cerr << "Use " << argv[0] << " -ffilename" << endl;
      exit(1);
    }
  while ((argc > 1) && (argv[1][0] == '-'))
    {
      switch (argv[1][1])
        {
        case 'f':
          filename = &argv[1][2];
          break;
        case 'm':
	  ss << argv[1][2];
	  ss >> m_id;
          break;
        case 'p':
	  print_out = true;
          break;
        default:
          cerr << "Bad option " << argv[1] << endl;
          exit(1);
        }
      ++argv;
      --argc;
    }

  return 0;
}
