//Get reference of the polyPatch at the inlet 
const polyPatch& pPatch =  mesh.boundaryMesh()[patchName];

const PrimitivePatch<faceList, pointField> priPatch
(
    pPatch.localFaces(), 
    pPatch.localPoints()
); 

//The big mesh file is not needed anymore 
mesh.clearOut();


const scalar mergeDist = 1.0e-10;
pointField points; 
faceList faces; 
labelList pointMap;

//Merge faces and points in the master
Foam::PatchTools::gatherAndMerge
(
    mergeDist,
    priPatch,
    points,
    faces,
    pointMap
);

pointField faceCentres; 
scalarField faceAreas; 

if(Pstream::master())
{
   faceCentres = pointField(faces.size());
   faceAreas = scalarField(faces.size());

   forAll(faces, facei)
   {
      faceCentres[facei] = faces[facei].centre(points);
      faceAreas[facei] = faces[facei].mag(points);
   }
}

Pstream::scatter(faceCentres);
Pstream::scatter(faceAreas);

Pout<< sum(faceAreas) << nl;


// // //Read the face centers  and area at the inlet on each processor 
// // // and collect them to the master processor.
// // //Use the patch name to extract from the boundary mesh
// // List<vectorField> gatheredPoints(Pstream::nProcs());
// // List<scalarField> gatheredAreas(Pstream::nProcs());

// // gatheredPoints[Pstream::myProcNo()] = mesh.boundary()[patchName].Cf();
// // gatheredAreas[Pstream::myProcNo()] = mesh.boundary()[patchName].magSf();

// // //Gather to the master 
// // Pstream::gatherList(gatheredPoints);
// // Pstream::gatherList(gatheredAreas);

// // //Scatter all the gathered to all the processors
// // Pstream::scatterList(gatheredPoints);
// // Pstream::scatterList(gatheredAreas);

// //Combine all to one list.
// vectorField points = ListListOps::combine<vectorField>(gatheredPoints, accessOp<vectorField>());
// scalarField areas = ListListOps::combine<scalarField>(gatheredAreas, accessOp<scalarField>());


// const label n = points.size();
// Info << "\nNumber of points on the patch = " << n << nl << nl;

// //Local points on the processor
// vectorField locPoints(0);

// //Local points on the processor
// labelList pointsIndx(0);

// //Fill the local points
// forAll(points, i)
// {
//    if(i%Pstream::nProcs() == Pstream::myProcNo())
//    {
//       locPoints.append(points[i]);
//       pointsIndx.append(i);
//    }
// }

// const label nLoc = locPoints.size();

// // Calculate the mean velocity, turbulence 
// // intensity, and length scale profiles. 
// vectorField Uav(nLoc, vector::zero);
// vectorField I(nLoc, vector::zero);
// vectorField L(nLoc, vector::zero);

// //Get Uav, I, L from the ABL profile.
// prof.getUav(locPoints, Uav);
// prof.getI(locPoints, I);
// prof.getL(locPoints, L);

// //Scale the turbulence intensity and length scales 
// //by the desired factors to account for decay of turbulence downstream.
// for(label i=0; i<nLoc; i++)
// {
//    I[i] = cmptMultiply(I[i], scaleI);
//    L[i] = cmptMultiply(L[i], scaleI);
// }

// //Frequency step
// scalar df = fMax/nFreq;	

// //The smallest frequency in the simulation
// scalar f0 = df/n;

// //Frequencies used for interpolating the CPSD matrix.
// scalarField fInterp(nInterp, 0.0);

// //The interpolation method is taken from Tao, Tianyou, et al. 
// // "Efficacy of Interpolation-Enhanced Schemes in Random Wind 
// //  Field Simulation over Long-Span Bridges." 
// //  Journal of Bridge Engineering 23.3 (2017): 0401714.
// forAll(fInterp,i)
// {
//      fInterp[i] = f0*::pow(n*nFreq, i/(nInterp - 1.0));
// }

// Info << "\nSelected interpolation frequencies = " << fInterp << nl << nl ;

// //Random number generator with a seed
// Random rand(seed);

// //Double indexed phase angles sign and consing computed.
// //Stored as a one dimensional array 
// vectorField sinPhi(nFreq*n, vector::zero);
// vectorField cosPhi(nFreq*n, vector::zero);


// scalar twopi = 2.0*Foam::constant::mathematical::pi;

// //Fill the double indexed phase angles
// if(Pstream::master())
// {      
//    //Pick  a vector random phase angels from a uniform distribution 
//    //in the interval [0.0, 2*pi]
//    for(int m = 0; m < n ; m++)
//    {
//       for(int l = 0; l < nFreq ; l++)
//       {
//          sinPhi[m*nFreq + l]  =  twopi*rand.sample01<vector>(); 
         
//          cosPhi[m*nFreq + l].x()  =  ::cos(sinPhi[m*nFreq + l].x()); 
//          cosPhi[m*nFreq + l].y()  =  ::cos(sinPhi[m*nFreq + l].y()); 
//          cosPhi[m*nFreq + l].z()  =  ::cos(sinPhi[m*nFreq + l].z());          
         
//          sinPhi[m*nFreq + l].x()  =  ::sin(sinPhi[m*nFreq + l].x()); 
//          sinPhi[m*nFreq + l].y()  =  ::sin(sinPhi[m*nFreq + l].y()); 
//          sinPhi[m*nFreq + l].z()  =  ::sin(sinPhi[m*nFreq + l].z()); 
//       }
//    }
// }

// Pstream::scatter(sinPhi);
// Pstream::scatter(cosPhi);

// //The coherency matrix is symmetric and only the lower triangular part is stored.
// //The matrix is stored in row-major packed Lapack format with n*(n+1)/2 entries 
// //including the diagonal elements
// scalarField Cu(n*(n+1)/2, 0.0);
// scalarField Cv(n*(n+1)/2, 0.0);
// scalarField Cw(n*(n+1)/2, 0.0);

// //List to store the decomposed the CPSD matrix for each interpolation frequency
// List<vectorField>H(n*nLoc, vectorField(nInterp, vector::zero));


// // mesh = NULL;



