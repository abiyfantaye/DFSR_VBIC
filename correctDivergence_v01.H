
const polyPatch infPatch= mesh.boundaryMesh()[patchName];

const pointField& patchPoints = infPatch.localPoints(); 

const List<face>& patchFaces  = infPatch.localFaces(); 


// const labelListList& edges = infPatch.faceEdges();

// SortableList<label> sortedPoints(patchPoints);

//All the points on new mesh back and front. 
pointField allPoints(2*patchPoints.size());

scalar dx = 0.1;

for(label i=0; i < patchPoints.size(); i++)
{
    allPoints[i] = patchPoints[i];
    allPoints[i + patchPoints.size()] = patchPoints[i] + vector(dx,0.0,0.0);
}

cellList allCells(0); 
faceList allFaces(0); 

Info << patchPoints << nl;
Info << patchFaces << nl;

label faceCount = 0;

forAll(patchFaces, faceI)
{
    //No. vertices for front and back faces. 
    label nVtxFront = patchFaces[faceI].size();

    //Number of faces = number of vertices + 2
    List<labelList> faces(nVtxFront + 2);

    faces[0] = labelList(nVtxFront);
    faces[1] = labelList(nVtxFront);
   
    //No. vertices for other faces. 
    label nVtxSides = 4;

    for(label i=0; i < nVtxFront; i++)
    {
        faces[0][i] = patchFaces[faceI][i];
        faces[1][i] = patchPoints.size() + patchFaces[faceI][i];
    }    

    for(label i=0; i < nVtxFront; i++)
    {
        faces[2+i] = labelList(nVtxSides);

        faces[2+i][0] = patchFaces[faceI][i];         
        faces[2+i][1] = patchFaces[faceI][i] +  patchPoints.size();      
        faces[2+i][2] = patchFaces[faceI][(i+1)%nVtxFront] +  patchPoints.size();        
        faces[2+i][3] = patchFaces[faceI][(i+1)%nVtxFront];        
    }
    
    labelList faceNames (faces.size(), -1);

    forAll(faces, i)
    {
        allFaces.append(face(faces[i]));    
        faceNames[i] = faceCount + i;
    }    
    
    faceCount += faces.size();

    allCells.append(cell(faceNames));

}


Info << allCells << nl;
Info << allFaces << nl;

fvMesh inletMesh
(
    Foam::IOobject
    (
        Foam::fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ
    ),
    pointField(allPoints),
    faceList(allFaces),
    cellList(allCells)
);

inletMesh.writeObject
( 
    IOstream::ASCII,
    IOstream::currentVersion,
    IOstream::UNCOMPRESSED
);


// const labelList meshPoints = mesh.boundary()[patchName].patch().meshPoints();
// const labelListList pointFaces = mesh.boundary()[patchName].patch().pointFaces();

// const labelListList faceFaces = mesh.boundary()[patchName].patch().faceFaces();

// const labelUList faceCells = mesh.boundary()[patchName].patch().faceCells();

// // Info << "Points = " <<  meshPoints << nl;
// // Info << "pointFace = " <<  pointFaces << nl;
// // Info << "Cells = " <<  faceCells << nl;
// Info << "faceFaces = " <<  faceFaces << nl;

// forAll(faceCells, i)
// {
//     Info << mesh.cellPoints(faceCells[i]) << nl;
// }


// mesh.boundary()[patchName].patch()


// Foam::fvMesh inletMesh
// (
//     Foam::IOobject
//     (
//         Foam::fvMesh::defaultRegion,
//         runTime.timeName(),
//         runTime,
//         Foam::IOobject::MUST_READ
//     )
// );