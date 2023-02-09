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

vectorField U(n, vector::zero);

scalar sqrt2df = ::sqrt(2.0*df);
scalar tolerance = 1.0e-10;
scalar f = 0.0;
vector Hjm;

//Generated velocity at each time step. 
for (int t = 0; t < Nt; t++)
{
     U = vector::zero;

     for (int r = 0; r < nLoc; r++)
     {
          label j = pointsIndx[r];
          for (int m = 0; m <= j; m++)
          {    
               for(int l = 0; l < nFreq; l++)
               {
                    f = (l + (m + 1.0)/n)*df;

                    Hjm = interpolateSpline(f, fInterp, H[r*n + m]);

                    if(mag(Hjm) > tolerance)
                    {                    
                         U[j].x() += abs(Hjm.x())*::cos(twopi*f*dt*t + Phi[l*n + m].x());
                         U[j].y() += abs(Hjm.y())*::cos(twopi*f*dt*t + Phi[l*n + m].y());
                         U[j].z() += abs(Hjm.z())*::cos(twopi*f*dt*t + Phi[l*n + m].z());
                    }    
               }
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
