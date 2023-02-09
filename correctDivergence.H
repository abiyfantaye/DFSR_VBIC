if(correctDivergence)
{   
    
    Ut -= flowDir*(sum(Ut.component(0)*faceAreas)/totArea);
    Ut += prof.getUav(faceCentres);

    label iteration  = 0;

    while(iteration < maxDivFreeIteration)
    {
        Info<< "DFSR: Calculating divergence corrections, Iteration = " << iteration << nl;

        scalar totSideFlux = 0.0;

        forAll(sidePatch, facei)
        {
            totSideFlux += sideAreas[facei]*(sideNormals[facei]&Ut[sideCells[facei]]);
        }
        
        //Set the net flux to zero to make it well posed  
        forAll(sidePatch, facei)
        {
            Ut[sideCells[facei]] -= (sideNormals[facei]*totSideFlux/totSideArea);
        }

        forAll(mergedFaces, celli)
        {
            U[celli] = Ut[celli];
        }    
        
        U.boundaryFieldRef()[0] = U.boundaryField()[0].patchInternalField();    
        
        lambdaBC = -U.boundaryField()[0] & sidePatch.faceNormals();
        
        Ux = U.component(0); 

        R = (fvc::ddt(Ux) + (fvc::grad(Ux)()&U))/Ux - fvc::div(U); //high order approx         
        // R = fvc::ddt(Ux)/Ux - fvc::div(U); //lower order approx    

        // R.primitiveFieldRef() = R*inletMesh.V();

        // fvc::smooth(R, 0.5);
        
        // forAll(mergedFaces, celli)
        // {
        //     R[celli] *= inletMesh.V()[celli]/dx;
        // }   

        R -= R.weightedAverage(inletMesh.V());        
        
        R.boundaryFieldRef()[0] = R.boundaryField()[0].patchInternalField();    

        fvScalarMatrix lambdaEqn 
        ( 
            fvm::laplacian(dimensionedScalar(dimless, 1), lambda) == R
        ); 

        lambdaEqn.solve();       

        volVectorField Ucorr("Ucorr", fvc::grad(lambda));        
        Ucorr.write();  

        U += fvc::grad(lambda);

        forAll(mergedFaces, celli)
        {
            Ut[celli] = U[celli];
        }    

        Info << "Continuity error: Initial = " << runTime.deltaTValue()*mag(R)().weightedAverage(inletMesh.V()).value()  
             << ", Final = " << runTime.deltaTValue()*mag(fvc::ddt(Ux)/Ux - fvc::div(U))().weightedAverage(inletMesh.V()).value() << nl;  
       
        iteration++; 
    }

    Ux.storePrevIter();

    inletMesh.write();        
    lambda.write();
    U.write();         
    R.write();          
}