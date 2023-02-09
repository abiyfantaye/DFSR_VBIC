Ut = vector::zero;

forAll(facesIndx, i)
{
    Ut[facesIndx[i]] = allU[runTime.timeIndex()][i];
}

reduce(Ut, sumOp<vectorField>());

if(correctDivergence)
{
    scalar totSideFlux = 0.0;

    forAll(sidePatch, facei)
    {
        totSideFlux += sideAreas[facei]*(sideNormals[facei]&Ut[sideCells[facei]]);
    }
    
    forAll(sidePatch, facei)
    {
       Ut[sideCells[facei]] -= (sideNormals[facei]*totSideFlux/totSideArea);
    }

    forAll(mergedFaces, celli)
    {
        U[celli] = Ut[celli]
    }


    // // Uf = fvc::interpolate(U); // & (inletMesh.Sf()/inletMesh.magSf());
    
    // // // vectorField patchU =  U.boundaryFieldRef()[0].patchInternalField();

    U.boundaryFieldRef()[0] = U.boundaryField()[0].patchInternalField();    
    
    lambdaBC = -U.boundaryField()[0] & sidePatch.faceNormals();


    // totSideFlux = 0.0;

    // forAll(sidePatch, facei)
    // {
    //     totSideFlux += sideAreas[facei]*(sideNormals[facei]&U.boundaryField()[0][facei]);
    // }    
    
    // Info << totSideFlux << nl;

    
    // forAll(sidePatch, facei)
    // {
    //     totFlux += sideAreas[facei]&U.boundaryField()[0][facei];
    // }

    // forAll(sidePatch, facei)
    // {
    //     totFlux += sideAreas[facei]*totFlux/totSideArea;
    // }



    // // vectorField patchU =  U.boundaryFieldRef()[0].patchInternalField();


    // // forAll(patch, facei)
    // // {
    // //     Info << patchU[facei]<< nl;
    // //     // Info << patch[facei] << nl;
    // //     // patchU[facei] = Uf.internalField()[patch.start() + facei];

    // //     // Info << Uf.internalField()[patch.start() + facei] << nl;
    // // }
    

    // // // fvc::interpolate(U)


    // // lambda = flowDir&U;

    // // volVectorField gradLambda1("gradLambda1", fvc::grad(lambda));
    
    // // gradLambda1.write();

    // // gradLambdaBC;

    // // forAll(mergedFaces, i)
    // // {
    // //     U[i] =  Ut[i];
    // // } 

    // // volScalarField R("R", flowDir&fvc::ddt(U));
    // // volScalarField divU("divU", fvc::div(U));

    // // volVectorField gradUx("gradUx", fvc::grad(vector(1,0,0)&U));

    // // volScalarField dUdt("dUdt", flowDir&fvc::ddt(U));
    // // volVectorField dUdx("dUdx", dimensionedScalar(dimless, 1.0)/(U + Uav));

    // // forAll(mergedFaces, i)
    // // {
    // //     Info << 1.0/(flowDir & U[i]) << nl;
    // // } 

    // // dUdx = (1.0/(flowDir&(Uav + U)))*dUdx;

    // // dUdx = 1.0/(flowDir&(Uav + U));

    // // volScalarField UavU("UavU", 1.0/(flowDir&(Uav + U)));

    // // volScalarField dUdx("dUdx", (flowDir&Uav));


    // // Info << min(Uav) << nl;
    // // Info << min(U) << nl;
    
    // // Uav.relax();


    // // Info << Uav.boundaryField() << nl;
    // // Info << U.boundaryField() << nl;

    // // volScalarField R1("R1", flowDir&fvc::ddt(U));
    
    // // volScalarField R2("R2", 1.0/);    
    volScalarField R("R", ((flowDir&fvc::ddt(U))/(flowDir&U))-fvc::div(U)); 

    // scalar totRFlux = 0.0;

    // forAll(faceAreas, facei)
    // {
    //     totRFlux += ;
    // }
    
    scalar sumAreaR = sum(faceAreas*f);
    
    forAll(faceAreas, facei)
    {
       R[facei] -= (sumAreaR/totArea);
       R[facei] = R[facei]*V[facei];
    }

    // Info << sum(faceAreas*R) << nl;


    // R.primitiveFieldRef() = faceAreas*R;

    // // R.primitiveFieldRef() = R.primitiveField() - sideFlux.primitiveField();

    // // scalar sumR = sum(R.primitiveField());

    // // R.primitiveFieldRef() -= sumR*faceAreas/totArea;


    // // Info << R << nl;
    // // Info << sum(R.primitiveField()) << nl;

    // // Info <<  lambda.boundaryField()[0] << nl;

    // // volScalarField R2("R2", 1.0/);    

    inletMesh.write();
    U.write();
    // // divU.write();
    // // gradUx.write();
    // // dUdx.write();
    // // R0.write();
    // // UavU.write();
    // // lambda.write();
    // U.write();    
    
    // // Uf.write();

    fvScalarMatrix lambdaEqn 
    ( 
        fvm::laplacian(dimensionedScalar(dimless, 1.0), lambda) == R 
    ); 

    // lambdaEqn.setReference(lambdaRefCell, lambdaRefValue);

    lambdaEqn.solve();    
    
    // lambda.correctBoundaryConditions();

    lambda.write();

    volVectorField gradLambda("gradLambda", fvc::grad(lambda));
    
    gradLambda.write();

    // volVectorField U0("U0", 1.0*U);

    // U0.write();

    // // // scalarField a1((U.boundaryField()[0] - gradLambda.boundaryField()[0])&(inletMesh.boundaryMesh()[0].faceNormals()));

    // // // Info << a1 << nl;
    // // // Info << inletMesh.boundaryMesh()[0].faceNormals() << nl;

    // Info<< "DFSR:  Applying divergence-free corrections" << nl;

    // if(Pstream::master())
    // {
    //     // Apply correction
    //     U+= fvc::grad(lambda);
    // }
    // // // volScalarField R("R", (flowDir&fvc::ddt(U))/(flowDir&U)-fvc::div(U));    

    // U.write();
    // R.write();

}

// //Add the mean to the fluctuations 
// if(Pstream::master())
// {
//     Ut += Uav.internalField();
// }