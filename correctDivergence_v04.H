if(correctDivergence)
{
    label iteration  = 0;

    while(iteration < maxDivFreeIteration)
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
            U[celli] = Ut[celli];
        }

        U.boundaryFieldRef()[0] = U.boundaryField()[0].patchInternalField();    
        
        lambdaBC = -U.boundaryField()[0] & sidePatch.faceNormals();
        
        Ux = U.component(0); 

        // phi = fvc::flux(U);
 
        // // u = U.component(0);  



        
        // if(Pstream::master())
        // {        
        //     volScalarField dudt("dudt", fvc::ddt(Ux));
        //     scalarField dudtMy = (1.5*u - u0 + 0.5*u00)/dt;
        //     // scalarField dudtMy = (u - u0)/dt;

        //     Info << "OF ddt = " << dudt[64] << nl;
        //     Info << "My ddt = " << dudtMy[64] << nl;
        // }
        
        // volScalarField dudx("dudx", fvc::ddt(Ux)/Ux);
        // volScalarField divU("divU", fvc::div(U));
        // volScalarField divPhi("divPhi", fvc::div(phi));

        // R = (fvc::ddt(Ux) + (fvc::grad(Ux)()&U))/Ux - fvc::div(U); //high order approx 
        R = fvc::ddt(Ux)/Ux - fvc::div(U); //low order approx 
        
        
        // forAll(mergedFaces, facei)
        // {
        //     R[facei] *= inletMesh.V()[facei];
        // }
        
        dimensionedScalar Rbar = R.weightedAverage(inletMesh.V());

        R -= Rbar;


        fvScalarMatrix lambdaEqn 
        ( 
            fvm::laplacian(dimensionedScalar(dimless, 1.0), lambda) == R
        ); 

        // lambdaEqn.setReference(-Pstream::myProcNo(), 0);
        
        // lambdaEqn.source() = R*inletMesh.V();    

        lambdaEqn.solve();    
       
        Info<< "DFSR: Applying divergence-free corrections" << nl;
        
        volVectorField Ucorr("Ucorr", fvc::grad(lambda));

        U += fvc::grad(lambda);

        // U.correctBoundaryConditions();

        forAll(mergedFaces, celli)
        {
            Ut[celli] = U[celli];
        }

        // Ux = U.component(0);  


        // volScalarField R1("R1", fvc::ddt(Ux)/Ux - fvc::div(U)); //low order approx 

        // volScalarField dudx1("dudx1", fvc::ddt(Ux)/Ux);
        // volScalarField divU1("divU1", fvc::div(U));


        // // R = (fvc::ddt(Ux) + (fvc::grad(Ux)()&U))/Ux - fvc::div(U); 
        // volScalarField R1("R1", fvc::ddt(Ux)/Ux - fvc::div(U)); //low order approx         
        
        Info << "Divergence error: " << runTime.deltaTValue()*mag(R)().weightedAverage(inletMesh.V()).value() << nl;           
        // Info << "Divergence error: " << runTime.deltaTValue()*mag(R1)().weightedAverage(inletMesh.V()).value() << nl;           
        

        inletMesh.write();        
        lambda.write();
        U.write();         
        // dudx.write(); 
        // divU.write();         
        // // divPhi.write();         
        // dudx1.write(); 
        // divU1.write(); 
        R.write(); 
        // R1.write(); 
        Ucorr.write();       
        iteration++; 
    }

    // u00 = u0;   
    // u0 = u;
}