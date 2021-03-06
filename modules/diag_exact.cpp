#include "diag_exact.hpp"

int diag_exact(int n, int nnz, complex<double> *nzval, int *irow, int *pcol)
{
  int cur_col, cur_row, ind_nzval, int_row;
  bool found = false;
  complex<double> *mat;
  complex<double> cur_el;
  gsl_vector_complex *eval;
  gsl_matrix_complex *evec;
  int i;

  mat = new complex<double>[n*n];

  ind_nzval = 0;
  for (cur_col = 0; cur_col < n; cur_col++)
     {
       for (cur_row = 0; cur_row < n; cur_row++)
 	{
	  found = false;
	  i = 0;
 	  for (int_row = irow[pcol[cur_col]] ;; int_row = irow[pcol[cur_col]+i])
 	    {
 	      if (int_row == cur_row)
 		{
		  if (ind_nzval < nnz)
		    insert_in_mat(mat, n*cur_row+cur_col, nzval[ind_nzval], n*n);
		  else
		    {
		      cerr << "Out of range index for nzval. " << endl;
		      exit(1);
		    }
 		  ind_nzval++;
 		  found = true;
		  break;
 		}
	      i++;
	      if (int_row == irow[pcol[cur_col+1]-1])
		break;
 	    }
	  if (found == false)
	    {
	      insert_in_mat(mat, n*cur_row+cur_col, 0, n*n);
	    }
	}
    }

  //cout << "Temporarily not printing matrix" << endl;
  cout << "****** Printing matrix ******* " << endl;
  cout << "{";
  for (cur_row = 0; cur_row < n; cur_row++)
    {
      cout << "{";
      for (cur_col = 0; cur_col < n; cur_col++)
  	{
	  cur_el = mat[n*cur_row+cur_col];
	  if (cur_col == n - 1)
	    {
	      if (cur_row == n - 1) //very last element
		{
		  if (imag(cur_el) == 0.)
		    cout << fixed << real(cur_el)  << "}} ";
		  else
		    cout << fixed << real(cur_el) << "+ I (" 
			 << fixed << imag(cur_el) << ") }} ";
		}
	      else
		{
		  if (imag(cur_el) == 0.)
		    cout << fixed << real(cur_el)  << "}, ";
		  else
		    cout << fixed << real(cur_el) << "+ I (" 
			 << fixed << imag(cur_el) << ") }, ";
		}
	    }
	  else
	    {
	      if (imag(cur_el) == 0.)
		cout << fixed << real(cur_el)  << ", ";
	      else
		cout << fixed << real(cur_el) << "+ I (" 
		     << fixed << imag(cur_el) << ") , ";
	    }
  	}
      cout << endl;
    }
  cout << "****************************** " << endl;

  /* Initialize eval, evec */
  eval = gsl_vector_complex_alloc (n);
  evec = gsl_matrix_complex_alloc (n, n);
  /* Diagonalize non symm matrix using gsl */
  diagonalize_non_symm(mat, n, eval, evec);

  /* Print result */
  cout << endl;
  cout << "Printing the eigenvalues " << endl;
  cout << endl;
  print_result(n, eval, evec);
  cout << endl;  

  /* Free */
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);
  
  delete [] mat;  

  return 0;
}

int insert_in_mat(complex<double> *mat, int pos, complex<double> val, int size)
{
  if (pos > size)
    {
      cerr << "Wrong assignement position in insert_in_mat. Exiting. " << endl;
      exit(1);
    }
  else
    {
      mat[pos] = val;
    }	

  return 0;
}

int diagonalize_non_symm(complex<double> *data, int size, gsl_vector_complex *eval,
                         gsl_matrix_complex *evec)
{
  int i;
  double data2[size*size];

  //convert from complex to double if elements are real
  for (i = 0; i < size*size; i++)
    {
      if (imag(data[i]) == 0)
	{
	  data2[i] = real(data[i]);
	}
      else
	{
	  cerr << "At the moment do not know what to do with complex non-hermitian." 
	       << "Exiting. " << endl;
	  exit(1);
	}
    }

  //If here it means that the matrix data2 is real and symmetric.

  gsl_matrix_view m
    = gsl_matrix_view_array (data2, size, size);

  gsl_eigen_nonsymmv_workspace * w =
    gsl_eigen_nonsymmv_alloc (size);

  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec,
                           GSL_EIGEN_SORT_ABS_DESC);

  return 0;
}

int print_result(int size, gsl_vector_complex *eval, gsl_matrix_complex *evec)
{ 
  int i;
  gsl_complex eval_i;
#ifdef SHOW_EIGENVECTORS
  gsl_vector_complex_view evec_i; 
  int j; 
#endif

  for(i = 0; i < size; i++)
    {
      eval_i = gsl_vector_complex_get (eval, i);

      cout << setprecision (10) << GSL_REAL(eval_i) << 
	" + " << setprecision (10) << GSL_IMAG(eval_i) << " I , " << endl;
#ifdef SHOW_EIGENVECTORS
      evec_i = gsl_matrix_complex_column (evec, i); 
      cout << "eigenvector = " << endl;
      for (j = 0; j < size; ++j)
	{ 
	  gsl_complex z = 
	    gsl_vector_complex_get(&evec_i.vector, j); 
	  cout <<  GSL_REAL(z) << " + " << GSL_IMAG(z) << " I " << endl; 
	} 
      cout << endl;
#endif
    }
  
  return 0;
}
