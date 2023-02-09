//- Returns the von Karman spectrum for i-th component given the
//  mean velocity(Uav), turbulence intensity(I) and length scales(L)
inline scalar vonK
(
  const label cmpt,
  const scalar f,
  const scalar Uav,
  const scalar I,
  const scalar L
)
{
  if(cmpt==0){
    return 4.0*::pow(I*Uav, 2.0)*(L/Uav)/::pow(1.0 + 70.8*::pow(f*L/Uav, 2.0), 5.0/6.0);
  }
  else if(cmpt==1){
    return 4.0*::pow(I*Uav, 2.0)*(L/Uav)*(1.0 + 188.4*::pow(2.0*f*L/Uav, 2.0))/::pow(1.0 + 70.8*::pow(2.0*f*L/Uav, 2.0), 11.0/6.0);
  }
  else if(cmpt==2){
    return 4.0*::pow(I*Uav, 2.0)*(L/Uav)*(1.0 + 188.4*::pow(2.0*f*L/Uav, 2.0))/::pow(1.0 + 70.8*::pow(2.0*f*L/Uav, 2.0), 11.0/6.0);
  }

  FatalErrorInFunction
  << "Component i = " << cmpt << " is unknown \n"
  << exit(FatalError);
  return 0;
}


//- Returns the coherency function for i-th component given freqquency,
//  coherency decay coefficient(C), mean velocity(Uav), and point pair(p1, p2).
inline scalar Coh
(  
    const label cmpt, 
    const scalar f,     
    const scalar Uav, 
    const tensor& C, 
    const vector& p1, 
    const vector& p2
)
{
  if(cmpt==0){
    return ::exp(-f*::sqrt(C.xx()*C.xx()*::pow(p2.x() - p1.x(), 2.0) + C.xy()*C.xy()*::pow(p2.y() - p1.y(), 2.0) + C.xz()*C.xz()*::pow(p2.z() - p1.z(), 2.0))/Uav);  
  }  
  else if(cmpt==1){
    return ::exp(-f*::sqrt(C.yx()*C.yx()*::pow(p2.x() - p1.x(), 2.0) + C.yy()*C.yy()*::pow(p2.y() - p1.y(), 2.0) + C.yz()*C.yz()*::pow(p2.z() - p1.z(), 2.0))/Uav); 
  }
  else if(cmpt==2){
    return ::exp(-f*::sqrt(C.zx()*C.zx()*::pow(p2.x() - p1.x(), 2.0) + C.zy()*C.zy()*::pow(p2.y() - p1.y(), 2.0) + C.zz()*C.zz()*::pow(p2.z() - p1.z(), 2.0))/Uav); 
  }
  
  FatalErrorInFunction
  << "Component i = " << cmpt << " is unknown \n"
  << exit(FatalError);
  return 0;
}

void printFullMatrix(int n, double *matrix)
{
  printf("\n \n");
  for (int i=0; i<n; i++)
  {
    for(int j=0; j<n; j++) 
    {
      printf("%.4f\t",matrix[i*n + j]);
    }
    printf("\n");
  }
}

void printFullMatrix(label n, label m, scalarField& matrix)
{
  printf("\n \n");
  for (int i=0; i<n; i++)
  {
    for(int j=0; j<m; j++) 
    {
      printf("%.4f\t",matrix[i*m + j]);
    }
    printf("\n");
  }
}

//Print the lower triangular part of a packed matrix(double)
void printPackedMatrix(label n, scalarField& matrix)
{
  printf("\n \n");
  for (int i=0; i<n; i++)
  {
    for(int j=0; j<=i; j++) 
    {
      printf("%.4f\t", matrix[j + i*(i+1)/2]);
    }
    printf("\n");
  }
}

//Check if the LDLT decomp has only 1 by 1 block diagonal decomposition
// -ipiv: is output from "LAPACKE_dsptrf" function that holds block information 
bool checkLDL(label n, labelField& ipiv)
{
  for (int i=0; i<n; i++)
  {
    if(ipiv[i] < 0)
    {
      return false; 
    }
  }
  return true;
}

//Computes the Cholesky factorization given LDL^T decomposition. 
//The lower triangular matrix of Cholesky decomposition is assigned to L.
void Cholesky(label n, scalarField& L)
{
  //List to hold the diagonal elements from LDL^T decomposition
  scalarField D(n, 0.0);

  //Fill the diagonals of the matrix.
  for (int i=0; i<n; i++)
  {
    //If Di is less than 0.0 the matrix is not PSD so 
    // compute the square root and L matrix only for posetive entries. 
    if(L[i + i*n] > 0.0)
    {
      D[i] = ::sqrt(L[i + i*n]);
    }
    else
    {
      D[i] = ::sqrt(L[(i-1) + (i-1)*n]);
    }
  }

  //Calculate the cholesky factorization 
  for (int i=0; i<n; i++)
  {
    L[i + i*n] = 1.0;

    for(int j=0; j<=i; j++)
    {
      L[j + i*n] *=  D[j];
    }
  }

  D.clear();
}

// //Computes the Cholesky factorization given LDL^T decomposition. 
// //The lower triangular matrix of Cholesky decomposition is assigned to L.
// void Cholesky(label n, scalarField& L)
// {
//   //List to hold the diagonal elements from LDL^T decomposition
//   scalarField D(n, 0.0);

//   //Fill the diagonals of the matrix.
//   for (int i=0; i<n; i++)
//   {
//     //If Di is less than 0.0 the matrix is not PSD so 
//     // compute the square root and L matrix only for posetive entries. 
//     if(L[i + i*(i+1)/2] > 0.0)
//     {
//       D[i] = ::sqrt(L[i + i*(i+1)/2]);
//     }
//   }

//   //Calculate the cholesky factorization 
//   for (int i=0; i<n; i++)
//   {
//     L[i + i*(i+1)/2] = 1.0;

//     for(int j=0; j<=i; j++)
//     {
//       L[j + i*(i+1)/2] *=  D[j];
//     }
//   }

//   D.clear();
// }

void fillLDL(label n, scalarField& d, scalarField& L)
{
  for( int i=0; i < n; i++)
  {
    d[i] = L[i + i*n];
    L[i + i*n] = 1.0;
  }  
}

// //Computes the Cholesky factorization given LDL^T decomposition. 
// //The lower triangular matrix of Cholesky decomposition is assigned to L.
// void Cholesky(label n, scalarField& L)
// {
//   scalar sum = 0.0;

//   for (label i = 0; i < n; ++i)
//   {
//     for (label j = 0; j <=i; ++j)
//     {
//       sum = L[j + i*(i+1)/2];

//       for (label k = 0; k < j; ++k)
//       {
//           sum -= L[k + i*(i+1)/2]*L[k + j*(j+1)/2];
//       }

//       if (i > j)
//       {
//           // L[j + i*(i+1)/2] = sum/L[j + j*(j+1)/2];
//       }
//       else if (sum > 0)
//       {
//           L[i + i*(i+1)/2] = ::sqrt(sum);
//       }
//       else
//       {
//         L[i + i*(i+1)/2] = 0.0;
//       }
//     }
//   }
// }