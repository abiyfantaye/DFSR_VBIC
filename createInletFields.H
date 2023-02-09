//Create fields at the inlet needed for div-free operation 
//Boundary conditions for lambda and U
wordList patchTypes(3); 

patchTypes[0] = "fixedGradient"; //sides dl/dx = 0, dl/dy = -v_c, dl/dz = -w_c
patchTypes[1] = "empty";  //front 
patchTypes[2] = "empty";  //back

//Lagrange multiplier to be solved
volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        inletMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    inletMesh,
    dimensionedScalar("lambda", dimensionSet(0,2,-1,0,0,0,0), 0.0), 
    patchTypes
); 

//Instantaneous velocity on all the cell/face centers 
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        inletMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    inletMesh,
    dimensionedVector("U", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
); 

//The source term.
volScalarField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        inletMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    inletMesh,
    dimensionedScalar("R", dimensionSet(0,0,-1,0,0,0,0), 0.0)
); 

// //Creating face fluxes.
// surfaceScalarField phi
// (
//     IOobject
//     (
//         "phi",
//         runTime.timeName(),
//         inletMesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     fvc::flux(U)
// ); 


//The velocity in x-direction.
volScalarField Ux
(
    IOobject
    (
        "Ux",
        runTime.timeName(),
        inletMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    U.component(0)
); 

// Ux.storePrevIter();


const polyPatch& sidePatch = inletMesh.boundaryMesh()["sides"];
const labelUList& sideCells = sidePatch.faceCells();
const scalarField sideAreas = sidePatch.magFaceAreas();
const vectorField sideNormals = sidePatch.faceAreas()/sideAreas;
const scalar totSideArea = sum(sideAreas);

//Boundary condition for lambda(fixedGradient)
fixedGradientFvPatchScalarField& lambadPatch = refCast<fixedGradientFvPatchScalarField>(lambda.boundaryFieldRef()[0]);
scalarField& lambdaBC =  lambadPatch.gradient();


