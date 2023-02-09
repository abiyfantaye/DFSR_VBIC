Info << "\nPerforming Cholesky decomposition for sampled frequencies" << nl << nl ; 

//Create the cross-power spectral density(CPSD) matrix. 

//Store the 2D coherency matrix into a column-major form for Cholesky factorization
//Since the matrix is symetric only the uper triangular part is stored
//The storage type is Lapack's 'packed' storage
float *Cu, *Cv, *Cw;

//Lists to store the lower triangular part of the decomposed matrix 
//List[frequency][row][column] format
List<List<vectorField>> H;

Cu = static_cast<float*>(mkl_malloc((n*(n+1)/2)*sizeof(float), 64));
Cv = static_cast<float*>(mkl_malloc((n*(n+1)/2)*sizeof(float), 64));
Cw = static_cast<float*>(mkl_malloc((n*(n+1)/2)*sizeof(float), 64));

if (Cu == NULL || Cv == NULL || Cw == NULL) {
  Info<< "\n ERROR: Can't allocate memory for matrices. Aborting... \n" << nl;
  mkl_free(Cu);
  mkl_free(Cv);
  mkl_free(Cw);
  return 1;
}

//spectra at each point.
vector s1(vector::zero);
vector s2(vector::zero);
vector coh(vector::zero);

//Tolerance 

scalar tol = 1.0e-6;

//Compute each entry of the matrix. 
//Since the matrix is symmetric only the upper triangular part is stored.
//Stored in row-major format
for( int f = 0; f < 1 ; f++)
{  
  for(int i=0; i < n; i++)
  {
    for (int j=0; j < n; j++)
    {
      s1 = vonKarmanSpectrum(fInterp[f], Uav[i], I[i], L[i]);
      s2 = vonKarmanSpectrum(fInterp[f], Uav[j], I[j], L[j]);
      
      coh = Coherency(fInterp[f], C, 0.5*(Uav[i] + Uav[j]), points[i], points[j]);

      // Su[n*i + j] = ::sqrt(s1.x()*s2.x())*coh.x();
      // Sv[n*i + j] = ::sqrt(s1.y()*s2.y())*coh.y();
      Cw[i + j*(j+1)/2] = coh.z();
      // Sw[n*i + j] = ::sqrt(s1.z()*s2.z())*coh.z();
      // Info << "i,j,k: " << i << "," << j << "," << i + j*(j+1)/2 << nl;

      // if(Su[n*i + j]< tol && i!=j)
      // {
      //   Su[n*i + j] = 0.0;
      // }
      // if(Sv[n*i + j]< tol && i!=j)
      // {
      //   Sv[n*i + j] = 0.0;
      // }
      // if(Sw[n*i + j]< tol && i!=j)
      // {
      //   Sw[n*i + j] = 0.0;
      // }
    }

    // //Increase the diagonal of the CPSD matrix by 
    // //some factor lambda to make it positive semi-defnite
    // Su[n*i + i]*=(lambda + 1.0);
    // Sv[n*i + i]*=(lambda + 1.0);
    // Sw[n*i + i]*=(lambda + 1.0);
    // Su[n*i + i]+=tol;
    // Sv[n*i + i]+=tol;
    // Sw[n*i + i]+=tol;
  }
  //Su *= 100.0;

// printPackedMatrix(n, Cw);


//  Info << "a(0,0) " << Sw[0] << nl;
//  Info << "a(0,1) " << Sw[1] << nl;
//  Info << "a(1,1) " << Sw[2] << nl;
//  Info << "a(0,2) " << Sw[3] << nl;
//  Info << "a(1,2) " << Sw[4] << nl;
//  Info << "a(2,2) " << Sw[5] << nl;

  // // // Perfrom Cholesky factorization using ScaLapack functions 
  // // int infoSu =  LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n , Su.data(), n);
  // // int infoSv =  LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n , Sv.data(), n);

  float *w1 = static_cast<float*>(mkl_malloc(n*sizeof(float), 64));
  float *w2 = static_cast<float*>(mkl_malloc(n*sizeof(float), 64));


  // int infoCw =  LAPACKE_spptrf(LAPACK_COL_MAJOR, 'U', n , Cw);
  
  mkl_sspffrt2(Cw, n, n,  'U', n , Cw);

  // if((infoCw != 0))
  // {
  //   FatalErrorInFunction
  //   << "Cholesky decomposition faild!\n"
  //   << "Please check you CPSD matrix or increase lambda.\n"
  //   << exit(FatalError);
  // }

  // if((infoSu != 0) || (infoSv != 0) || (infoSw != 0))
  // {
  //   FatalErrorInFunction
  //   << "Cholesky decomposition faild!\n"
  //   << "Please check you CPSD matrix or increase lambda.\n"
  //   << exit(FatalError);
  // }

  // //Copy the decomposed matrix into the global array
  // List<vectorField> Hf;
  // for(int i=0; i < n; i++)
  // {
  //   vectorField Hi(i + 1, vector::zero);

  //   for (int j=0; j <=i; j++)
  //   {
  //     Hi[j].x() = Su[n*i + j];
  //     Hi[j].y() = Sv[n*i + j];  
  //     Hi[j].z() = Sw[n*i + j]; 
  //   }
  //   Hf.append(Hi);
  // }

  // H.append(Hf);

  // printPackedMatrix(n, Cw);


  Info <<100.0*(f + 1)/nInterp << "\% completed" << ", f = " << fInterp[f] << nl << nl;
  

  // double* C = static_cast<double*>(mkl_calloc(n*n, sizeof(double), 64));
  // double* M = static_cast<double*>(mkl_calloc(n*n, sizeof(double), 64));

  // for(int i=0; i < n; i++)
  // {
  //   for (int j=i; j < n; j++)
  //   {
  //     C[i + j*n] = static_cast<double>(Cw[i + j*(j+1)/2]);
  //     //C[j*n + i] = C[i*n + j];
  //   }
  // }



  // // scalarField Sr(n*n, 0.0);

  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, C, n, C, n, 0.0, M, n);

  // printFullMatrix(n, n, M);

}
 
