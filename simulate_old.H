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
          pointsOs << points; 
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
                    if(mag(points[j] - samplingPoints[i]) < dist)
                    {
                         dist = mag(points[j] - samplingPoints[i]);
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
vectorField U(n, vector::zero);

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

for (int t = 0; t < Nt; t++)
{
     U = vector::zero;

     for (int r = 0; r < nLoc; r++)
     {
          label j = pointsIndx[r];
          for (int m = 0; m <= j; m++)
          {    
               for(int l = 0; l < M; l++)
               {
                    if(l < nFreq)
                    { 
                         const vector& Hjm = interpolateSpline((l + (m + 1.0)/n)*df, fInterp, H[r*n + m]);

                         Bu[l][0] = Hjm.x()*::cos(Phi[l*n + m].x());
                         Bv[l][0] = Hjm.y()*::cos(Phi[l*n + m].y());    
                         Bw[l][0] = Hjm.z()*::cos(Phi[l*n + m].z());
                         

                         Bu[l][1] = Hjm.x()*::sin(Phi[l*n + m].x());
                         Bv[l][1] = Hjm.y()*::sin(Phi[l*n + m].y());
                         Bw[l][1] = Hjm.z()*::sin(Phi[l*n + m].z());
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

               scalar temp = t*dt*2*pi*(m + 1.0)*df/n;
               scalar c = ::cos(temp);
               scalar s = ::sin(temp);
               label q = t%M;

               U[j].x() += (Du[q][0]*c - Du[q][1]*s); 
               U[j].y() += (Dv[q][0]*c - Dv[q][1]*s); 
               U[j].z() += (Dw[q][0]*c - Dw[q][1]*s); 
          }    

          //Multiply by sqrt(2*df)
          U[j]*=sqrt2df;   
     }
     

     
     //Add the mean. 
     forAll(pointsIndx, i)
     {
          U[pointsIndx[i]] += Uav[i];
     }

     reduce(U, sumOp<vectorField>());

     if(Pstream::master())
     {    
          if(writeInflow)
          {
               fileName timePath = inflowPath/name(t*dt);
               mkDir(timePath);

               OFstream uOS(timePath/"U");
               uOS << U;
          }

          if(writeSamples)
          {                
               sampleUxOS.open(samplePath/"Ux", std::ios_base::app);
               sampleUyOS.open(samplePath/"Uy", std::ios_base::app);
               sampleUzOS.open(samplePath/"Uz", std::ios_base::app);   

               sampleUxOS <<  t*dt << "\t";
               sampleUyOS <<  t*dt << "\t";
               sampleUzOS <<  t*dt << "\t";

               forAll(samplingPoints, i)
               {
                    sampleUxOS << U[samplingIndx[i]].x() << "\t";
                    sampleUyOS << U[samplingIndx[i]].y() << "\t";
                    sampleUzOS << U[samplingIndx[i]].z() << "\t";
               }

               sampleUxOS << nl;
               sampleUyOS << nl;
               sampleUzOS << nl;  

               sampleUxOS.close();
               sampleUyOS.close();
               sampleUzOS.close();
          }
     }     
     Info<< "Time = " << t*dt << nl;
     Info << 100.0*(t + 1.0)/Nt << "\% completed" << nl
          << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
          << "  ClockTime = " << runTime.elapsedClockTime() << " s"
          << nl << endl;
}

fftw_free(Bu); fftw_free(Du); 
fftw_free(Bv); fftw_free(Dv); 
fftw_free(Bw); fftw_free(Dw);
