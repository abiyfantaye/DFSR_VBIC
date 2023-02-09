Info << "\nPerforming Cholesky decomposition for sampled frequencies" << nl << nl ; 

// Info <<"Created CPSD matrix" << nl << nl;

//A List to store the CPSD matrix for each interplation frequency
//Since the matrix is symmetric only the lower triangular part is stored.
//The CPSD matrix is stored in row-major packed Lapack format. 
//The lower triangular part of the matrix has n*(n+1)/2 entries 
//including the diagonal of the matrix
List<scalarField>cpsd(nInterp, scalarField(n*(n+1)/2, 0.0));

scalarField cpsdCopy(n*(n+1)/2, 0.0);

//The spectra at a particular frequency
scalarField S(n, 0.0);


label cmpt = 2;

//The Cholesky decomposition is done using the following simplified form 
/*
  CPSD = diag(S)*C*diag(S)
  CPSD = diag(S)*Lc*Lc^T*diag(S)

  Where: 
        Lc: is the Cholesky factorization of the coherency matrix
        S: is the auto-power spectrum at each point(size = n)
        diag(S) is a diagonal matrix created from S. 
  
  Thus, the Cholesky factorization of the CPSD matrix is computed as:

  L =  diag(S)*Lc

  Note: This approach prevents the Cholesky factorization from failing due to 
        non-positive semi-defnitness of the CPSD matrix caused by the inhomogeneity 
        in the coherency matrix(covariance matrix) 
*/

//Compute each entry of the coherency matrix.
forAll(fInterp, f)
{  
  forAll(points, i)
  {
    S[i] = vonK(cmpt, fInterp[f], Uav[i], I[i][cmpt], L[i][cmpt]);

    for (int j=0; j <= i; j++)
    {
      cpsd[f][j + i*(i+1)/2] = Coh(cmpt, fInterp[f], cohUav, C, points[i], points[j]); 
      // cpsdCopy[j + i*(i+1)/2] = Coh(cmpt, fInterp[f], cohUav, C, points[i], points[j]); 

      cpsdCopy[j + i*(i+1)/2] = ::sqrt(vonK(cmpt, fInterp[f], Uav[i], I[i][cmpt], L[i][cmpt])
                                      *vonK(cmpt, fInterp[f], Uav[j], I[j][cmpt], L[j][cmpt]))
                                      *Coh(cmpt, fInterp[f], cohUav, C, points[i], points[j]);
    }  
  }

  //Compute the Cholesky factorization using LPACKE subroutine 
  //Uses a packed matrix storage.
  int info =  LAPACKE_dpptrf(LAPACK_ROW_MAJOR, 'L', n , cpsd[f].data());

  //If decomposition fail, show error message.
  if((info != 0))
  {
    FatalErrorInFunction
    << "ERROR: Cholesky decomposition failed for component = " << cmpt << nl
    << exit(FatalError);
  }

  CholeskyCPSD(n, S, cpsd[f]);

  Info  << 100.0*(f + 1)/nInterp << "\% completed" << ", f = " << fInterp[f] << nl
        << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;
  
  // printFullMatrix(n, n, cpsd[f]);


  scalarField C(n*n, 0.0);
  scalarField M(n*n, 0.0);
  

  printPackedMatrix(n, cpsdCopy);

  for(int i=0; i < n; i++)
  {
    for (int j=0; j<=i; j++)
    {
      C[j + i*n] = cpsd[f][j + i*(i+1)/2];
    }
  }


  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1.0, C.data(), n, C.data(), n, 0.0, M.data(), n);  
    


  // M -= C;

  printFullMatrix(n, M);

  // for(int i=0; i < n; i++)
  // {
  //   // Info << M[i + i*n] << nl;
  //   printf("%.6f\t", M[i + i*n]);
  //   printf("%.6f\t", cpsdCopy[i + i*n]);
  //   printf("%.6f\t", cpsd[f][i + i*n]);
  //   Info << nl;
  // }

}
 
