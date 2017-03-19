#include "diag_arn.hpp"

int diag_arn(int n, int nnz, double *nzval, int *irow, int *pcol, int neval,
	     bool print_out, std::ofstream &outfile)
{
  ARluNonSymMatrix<double, double> matrix(n, nnz, nzval, irow, pcol);

  // Defining what we need: neval eigenvectors matrix A with largest
  // magnitude (LM)
  ARluNonSymStdEig<double> dprob(neval, matrix);
  dprob.ChangeWhich("SR");
  //dprob.ChangeWhich("LR"); // LARGEST REAL
  // Turns on the debug mode
  dprob.Trace();
  // dprob.ChangeNcv(2*dprob.GetNev()+1);
  // Finding eigenvalues and eigenvectors.
#ifdef SHOW_EIGENVECTORS
  dprob.FindEigenvectors();
#endif
#ifndef SHOW_EIGENVECTORS
  dprob.FindEigenvalues();
#endif

  // Printing solution.
  Solution(matrix, dprob);

  if (print_out == true)
    {
      //Printing eigenvectors on outfile
      Eigenv_outfile(matrix, dprob, outfile);
    }

  return 0;
}
