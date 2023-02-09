     
     
     const fileName procPath(args.path()/"constant"/"SRFGData");


     mkDir(procPath);

     //Writes the inlet velocity field to a file. 
     OFstream UOs (procPath/"U", IOstream::BINARY);

     UOs << Ugen;
    
     //Pout <<"Processor "<< Pstream::myProcNo() << " finished writing!" << nl << nl; 
