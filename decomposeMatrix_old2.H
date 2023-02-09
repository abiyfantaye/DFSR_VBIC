Info << "\nPerforming Cholesky decomposition for sampled frequencies" << nl << nl ; 

// Info <<"Created CPSD matrix" << nl << nl;

//Lists to store CPSD matrix for each interplation frequency
//Since the matrix is symmetric only the lower triangular part is stored.
//The matrix is stored in row-major format packed Lapack format
List<scalarField>cpsd(nInterp, scalarField(n*(n+1), 0.0));

label cmpt = 2;

//Compute each entry of the matrix. 
for(int f = 0; f < nInterp; f++)
{  
  for(int i=0; i < n; i++)
  {
    for (int j=0; j <= i; j++)
    {
      // cpsd[f][j + i*(i+1)/2] = Coh(cmpt, fInterp[f], 0.5*(Uav[i] + Uav[j]), C, points[i], points[j]);
      cpsd[f][j + i*(i+1)/2] = Coh(cmpt, fInterp[f], cohUav, C, points[i], points[j]);

      // cpsd[f][j + i*(i+1)/2] = ::sqrt(vonK(cmpt, fInterp[f], Uav[i], I[i][cmpt], L[i][cmpt])
      //                                *vonK(cmpt, fInterp[f], Uav[j], I[j][cmpt], L[j][cmpt]))
      //                                 *Coh(cmpt, fInterp[f], 15.0, C, points[i], points[j]);
    }
    
    // cpsd[f][i + i*n]*(1)
    
  }  
  
  // int info =  LAPACKE_dsptrf(LAPACK_ROW_MAJOR, 'L', n , cpsd[f].data(), ipiv.data());

    // mkl_sspffrt2

  diag()  
  scalarField cpsdCopy(cpsd[f]);
  // for(int i=0; i < n; i++)
  // {
  //   for (int j=0; j<=i; j++)
  //   {
  //     cpsdCopy[j + i*n] = cpsd[f][j + i*n];
  //   }
  // }

  // printPackedMatrix(n, cpsd[f]);

  int info =  LAPACKE_dpotrf(LAPACK_ROW_MAJOR, 'L', n , cpsd[f].data(), n);


  // Info << info << nl;

  // int info =  LAPACKE_dpptrf(LAPACK_ROW_MAJOR, 'L', n , cpsd[f].data());

  // if(info != 0)
  // {

    // Declear working arrays 
    labelField ipiv(n, 0);
    scalarField e(n, 0.0);
    scalarField d(n, 0.0);
    
    // Perform LDL decomposition with for symetric indefnite matrix
    // int infoChol =  LAPACKE_dsytrf(LAPACK_ROW_MAJOR, 'L', n , cpsd[f].data(), n, ipiv.data());

    // int infoConv =  LAPACKE_dsyconv(LAPACK_ROW_MAJOR, 'L', 'R', n , cpsd[f].data(), n, ipiv.data(), e.data());

    // infoConv =  LAPACKE_dsyconv(LAPACK_ROW_MAJOR, 'L', 'R', n , cpsd[f].data(), n, ipiv.data(), e.data());
    
    // Cholesky(n, cpsd[f]);

    //Check if the decomposition succeded 
    if((info != 0))
    {
      FatalErrorInFunction
      << "ERROR: Cholesky decomposition faild for component = " << cmpt << nl
      << exit(FatalError);
    }




    // //Check if the LDL^T decomposition has the correct form(1x1 block-diagonal)
    // if(!checkLDL(n,ipiv))
    // {
    //   FatalErrorInFunction
    //   << "LDL^T decomposition has incompatable format!\n"
    //   << exit(FatalError);
    // }

  //   Cholesky(n, cpsd[f]);
  // }


  // fillLDL(n, d, cpsd[f]);



  // int info = LAPACKE_dsterf(n, d.data(), e.data());

  // //Check if the decomposition succeded 
  // if(info != 0)
  // {
  //   Info << info << nl;
  //   FatalErrorInFunction
  //   << "ERROR: Cholesky decomposition faild for component = " << cmpt << nl
  //   << exit(FatalError);
  // }  
  
  // for(int i=0; i < n; i++)
  // {
  //   Info << e[i] << nl; 
  //   // if (d[i] < 0.0)
  //   // {
  //   //   d[i] = 0.0;
  //   // }

  //   // if (e[i] < 0.0)
  //   // {
  //   //   d[i] = 0.0;
  //   // }
  // }

  // for(int i=0; i < n; i++)
  // {
  //   Info << d[i] << ", "; 
  //   // if(i < n-1)
  //   Info << e[i]<< nl; 
  //   // printf("%.10f\n", d[i]);

  // }




  Info <<100.0*(f + 1)/nInterp << "\% completed" << ", f = " << fInterp[f] << nl << nl;
  
  // printFullMatrix(n, n, cpsd[f]);


  // double* C = static_cast<double*>(mkl_calloc(n*n, sizeof(double), 64));
  // double* M = static_cast<double*>(mkl_calloc(n*n, sizeof(double), 64));

  // for(int i=0; i < n; i++)
  // {
  //   for (int j=0; j<=i; j++)
  //   {
  //     C[j + i*n] = cpsd[f][j + i*n];
  //     //C[j*n + i] = C[i*n + j];
  //   }
  // }


  // // int infoLLT = LAPACKE_dlauum(LAPACK_ROW_MAJOR, 'L', n, C, n);

  // // Info << infoLLT << nl;

  // printFullMatrix(n, C);


  // scalarField C(n*n, 0.0);
  // scalarField M(n*n, 0.0);
  

  // printFullMatrix(n, n, cpsdCopy);

  // for(int i=0; i < n; i++)
  // {
  //   for (int j=0; j<=i; j++)
  //   {
  //     C[j + i*n] = cpsd[f][j + i*(i+1)/2];
  //   }
  // }


  // cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, cpsd[f].data(), n, cpsd[f].data(), n, 0.0, M.data(), n);  
    
  // int infoLLT = LAPACKE_dlauum(LAPACK_ROW_MAJOR, 'L', n, cpsd[f].data(), n);


  // M = M - cpsdCopy;

  // printFullMatrix(n, n, M);

  // for(int i=0; i < n; i++)
  // {
  //   // Info << M[i + i*n] << nl;
  //   printf("%.6f\t", M[i + i*n]);
  //   printf("%.6f\t", cpsdCopy[i + i*n]);
  //   printf("%.6f\t", cpsd[f][i + i*n]);
  //   Info << nl;
  // }

}
 
