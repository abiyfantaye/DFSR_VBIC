Info <<"\n\nSimulating velocity time history:  " << nl << nl ;

fileName inflowPath = runTime.constant()/"boundaryData"/patchName;
fileName samplePath = runTime.constant()/"boundaryData"/"windProfile"/"sampledData";

vectorField samplingPoints;
labelList samplingIndx;

ofstream sampleUxOS;
ofstream sampleUyOS;
ofstream sampleUzOS;

//Create the directory where the inflow 
//will be stored. 
if(Pstream::master())
{
     if(writeInflow)
     {
          mkDir(inflowPath);

          OFstream pointsOs(inflowPath/"points");
          pointsOs << faceCentres; 
     }

     if(writeSamples)
     {
          mkDir(samplePath);
          IFstream sampleIS(samplePath/"samplingPoints");
          sampleIS >> samplingPoints; 

          samplingIndx = labelList(samplingPoints.size(), -1);
          scalar dist = VGREAT;

          //Find the index of the nearest points to the samplingPoints
          for(label i = 0; i < samplingPoints.size(); i++)
          {
               for(int j = 0; j < n ; j++)
               {
                    if(mag(faceCentres[j] - samplingPoints[i]) < dist)
                    {
                         dist = mag(faceCentres[j] - samplingPoints[i]);
                         samplingIndx[i] = j;
                    }
               }
               dist = VGREAT;
          }
          
          sampleUxOS.open(samplePath/"Ux");
          sampleUyOS.open(samplePath/"Uy");
          sampleUzOS.open(samplePath/"Uz");

          // //Write the nearest points 
          if(sampleUxOS.good() && sampleUyOS.good() && sampleUzOS.good())
          {
               sampleUxOS << "# Nearest points (x, y, z)" << nl;
               sampleUyOS << "# Nearest points (x, y, z)" << nl;
               sampleUzOS << "# Nearest points (x, y, z)" << nl;          
               
               forAll(samplingPoints, i)
               {
                    sampleUxOS << "# Point " << i << " (" << points[samplingIndx[i]].x() << ", " << points[samplingIndx[i]].y() << ", " << points[samplingIndx[i]].z() << ")" << nl;
                    sampleUyOS << "# Point " << i << " (" << points[samplingIndx[i]].x() << ", " << points[samplingIndx[i]].y() << ", " << points[samplingIndx[i]].z() << ")" << nl;
                    sampleUzOS << "# Point " << i << " (" << points[samplingIndx[i]].x() << ", " << points[samplingIndx[i]].y() << ", " << points[samplingIndx[i]].z() << ")" << nl;
               }
          }

          sampleUxOS.close();
          sampleUyOS.close();
          sampleUzOS.close();
     }

}

//Generated velocity at the current time. 
List<vectorField> U(Nt, vectorField(nLoc, vector::zero));

//Number of frequencies in the FFT 
label M(2*nFreq);

//Arrays to hold inverse FFT variables.
fftw_complex *Bu, *Bv, *Bw, *Du, *Dv, *Dw;


Bu = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));
Du = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));

Bv = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));
Dv = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));

Bw = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));
Dw = static_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex) * M));

fftw_plan pu;
fftw_plan pv;
fftw_plan pw;

scalar sqrt2df = ::sqrt(2.0*df);
vector Hjm;
scalar f,temp, cs, sn;
label q;

for (label i = 0; i < nLoc; i++)
{
     label j = pointsIndx[i];

     for (label m = 0; m <= j; m++)
     {    
          for(label l = 0; l < M; l++)
          {
               if(l < nFreq)
               { 
                    Hjm = sqrt2df*interpolateSpline((l + (m + 1.0)/n)*df, fInterp, H[i*n + m]);

                    Bu[l][0] = Hjm.x()*cosPhi[l*n + m].x();
                    Bv[l][0] = Hjm.y()*cosPhi[l*n + m].y();    
                    Bw[l][0] = Hjm.z()*cosPhi[l*n + m].z();
                    

                    Bu[l][1] = Hjm.x()*sinPhi[l*n + m].x();
                    Bv[l][1] = Hjm.y()*sinPhi[l*n + m].y();
                    Bw[l][1] = Hjm.z()*sinPhi[l*n + m].z();
               }
               else
               {
                    Bu[l][0] = 0.0; Bu[l][1] = 0.0;
                    Bv[l][0] = 0.0; Bv[l][1] = 0.0;
                    Bw[l][0] = 0.0; Bw[l][1] = 0.0;
               }
          }

          pu = fftw_plan_dft_1d(M, Bu, Du, FFTW_BACKWARD, FFTW_ESTIMATE);
          pv = fftw_plan_dft_1d(M, Bv, Dv, FFTW_BACKWARD, FFTW_ESTIMATE);
          pw = fftw_plan_dft_1d(M, Bw, Dw, FFTW_BACKWARD, FFTW_ESTIMATE);

          fftw_execute(pu);
          fftw_execute(pv);
          fftw_execute(pw);

          fftw_destroy_plan(pu);
          fftw_destroy_plan(pv);
          fftw_destroy_plan(pw);

          f = dt*twopi*(m + 1.0)*df/n;

          for (label t = 0; t < Nt; t++)
          {
               temp = t*f;
               cs = ::cos(temp);
               sn = ::sin(temp);
               q = t%M;

               U[t][i].x() += (Du[q][0]*cs - Du[q][1]*sn); 
               U[t][i].y() += (Dv[q][0]*cs - Dv[q][1]*sn); 
               U[t][i].z() += (Dw[q][0]*cs - Dw[q][1]*sn); 
          }  
     } 

     Info << "Generation completed for "<< i*Pstream::nProcs() << " points out of " << n << nl
          << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << nl << endl;
    
}

cosPhi.clear();
sinPhi.clear();
H.clear();

//Correct mass balance at the inlet
if(correctDivergence)
{
     scalar totArea = sum(faceAreas);
     vector flowDir(1.0, 0.0, 0.0);
     scalar massFlux;

     for(label t = 0; t < Nt; t++)
     {    
          massFlux = 0.0;

          forAll(facesIndx, i)
          {
               massFlux += (flowDir & U[t][i])*faceAreas[facesIndx[i]];
          }

          reduce(massFlux, sumOp<scalar>());

          //Correct to have sum(Ai*Uxi) = 0. 
          U[t] -= vector(massFlux/totArea, 0.0, 0.0);
     }     
     
     Info <<"\n\nApplyied mass flux correction" << nl << nl ;
}

vectorField gUt(n, vector::zero);

//Add to the global array. 
for(label t = 0; t < Nt; t++)
{     
     //Add the mean to the fluctuation. 
     U[t] += Uav;

     gUt = vector::zero;

     forAll(pointsIndx, i)
     {
          gUt[pointsIndx[i]] =  U[t][i];
     }

     reduce(gUt, sumOp<vectorField>());

     if(Pstream::master())
     {    
          if(writeInflow)
          {
               fileName timePath = inflowPath/name(startTime + t*dt);
               mkDir(timePath);

               OFstream uOS(timePath/"U");
               uOS << gUt;
          }

          if(writeSamples)
          {                
               sampleUxOS.open(samplePath/"Ux", std::ios_base::app);
               sampleUyOS.open(samplePath/"Uy", std::ios_base::app);
               sampleUzOS.open(samplePath/"Uz", std::ios_base::app);   

               sampleUxOS <<  startTime + t*dt << "\t";
               sampleUyOS <<  startTime + t*dt << "\t";
               sampleUzOS <<  startTime + t*dt << "\t";

               forAll(samplingPoints, i)
               {
                    sampleUxOS << gUt[samplingIndx[i]].x() << "\t";
                    sampleUyOS << gUt[samplingIndx[i]].y() << "\t";
                    sampleUzOS << gUt[samplingIndx[i]].z() << "\t";
               }

               sampleUxOS << nl;
               sampleUyOS << nl;
               sampleUzOS << nl;  

               sampleUxOS.close();
               sampleUyOS.close();
               sampleUzOS.close();
          }
     }

     Info<< "Writing time = " << startTime + t*dt << nl;

}

fftw_free(Bu); fftw_free(Du); 
fftw_free(Bv); fftw_free(Dv); 
fftw_free(Bw); fftw_free(Dw);
