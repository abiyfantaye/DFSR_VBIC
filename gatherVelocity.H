//Instantaneous velocity at time t
//gather from all processors 
vectorField Ut(n, vector::zero);

forAll(facesIndx, i)
{
    Ut[facesIndx[i]] = allU[runTime.timeIndex()][i];
}

reduce(Ut, sumOp<vectorField>());
   
