const PrimitivePatch<faceList, pointField> infPatch
(
    mergedFaces, 
    mergedPoints
); 

const IOobject io
(
    Foam::fvMesh::defaultRegion,
    runTime.timeName(),
    runTime,
    Foam::IOobject::MUST_READ
);

extrudedPatch infExtrude(io, infPatch);

fvMesh infFvMesh
(
    io,
    pointField(infExtrude.points()),
    faceList(infExtrude.faces()),
    cellList(infExtrude.cells())
);

infFvMesh.addFvPatches(infExtrude.patches(infPatch));

// Info << infFvMesh.boundaryMesh() << nl;
// infFvMesh.write();


// vectorField U(n, vector::zero);

//Boundary conditions for lambda
wordList patchTypes(3); 
patchTypes[0] = "zeroGradient";
patchTypes[1] = "empty";
patchTypes[2] = "empty"; 

volScalarField lambda
(
    IOobject
    (
        "lambda",
        runTime.timeName(),
        infFvMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    infFvMesh,
    dimensionedScalar("lambda", dimensionSet(0,2,-1,0,0,0,0), 0.0), 
    patchTypes
); 

volVectorField Uyz
(
    IOobject
    (
        "Uyz",
        runTime.timeName(),
        infFvMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    infFvMesh,
    dimensionedVector("Uyz", dimensionSet(0,1,-1,0,0,0,0), vector::zero)
); 

volScalarField Ux
(
    IOobject
    (
        "Ux",
        runTime.timeName(),
        infFvMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    infFvMesh,
    dimensionedScalar("Ux", dimensionSet(0,1,-1,0,0,0,0), 0.0)
); 

volScalarField Uxav
(
    IOobject
    (
        "Uxav",
        runTime.timeName(),
        infFvMesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    infFvMesh,
    dimensionedScalar("Uxav", dimensionSet(0,1,-1,0,0,0,0), 0.0)
); 

Info<< "\nApplying divergence-free corrections\n" << endl;

while (runTime.loop())
{    
    Info<< "Time = " << runTime.timeName() << nl << endl;

    forAll(mergedFaces, i)
    {
        Uyz[i] = U[runTime.timeIndex()][i];
        Uyz[i].x() = 0.0;
        Ux[i] = U[runTime.timeIndex()][i].x();
    }

    // R = fvc::ddt(Ux)/(Uxav + Ux) - fvc::div(Uyz);
    fvScalarMatrix lambdaEqn 
    ( 
        fvm::laplacian(dimensionedScalar(dimless, 1.0), lambda) 
        == 
        fvc::ddt(Ux)/(Uxav + Ux) - fvc::div(Uyz)
    ); 

    lambdaEqn.solve();

    Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

}
// for(label t=0; t<Nt; t++)
// {


//     //Use backward differencing scheme for time derivaties 
//     // dudt = (1.0/dt)*(1.5*Ux - 2.0*Ux0 + 0.5*Ux00);     
    
//     //First order approach 
//     dudt = (Ux - Ux0)/dt; 


//     





//     Info<< R << nl;

// }




// // infFvMesh.write();

// volVectorField volU
// (
//     IOobject
//     (
//         "volU",
//         runTime.timeName(),
//         infFvMesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     infFvMesh,
//     dimensionedVector("R", dimensionSet(0,0,-1,0,0,0,0), vector::zero)
// ); 

// volScalarField volU(U);

// faceCentres divU
// (
//     IOobject
//     (
//         "divU",
//         runTime.timeName(),
//         infFvMesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     infFvMesh,
//     dimensionedScalar("divU", dimensionSet(0,0,-1,0,0,0,0), 0.0), 
//     patchTypes
// ); 


//Assign boundary conditions 


// Info << lambda.boundaryField().types() << nl;
// Info << lambda << nl;




// fvOptions.constrain(TEqn); 
// TEqn.solve(); 
// fvOptions.correct(T); 

// Info<< "End\n" << endl;





// infPolyMesh.clearOut();

// infFvMesh.write();

