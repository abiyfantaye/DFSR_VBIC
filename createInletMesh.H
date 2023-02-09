
//Reverse face normals so that the extruded mesh 
//extends in x-direction 
forAll(mergedFaces, facei)
{
    mergedFaces[facei] = reverseList(mergedFaces[facei]);
}

const PrimitivePatch<faceList, pointField> inletPatch
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

//Approximate mesh size in x-direction 
const scalar dx = Foam::sqrt(totArea/faceCentres.size());

extrudedPatch infExtrude(io, inletPatch, dx);

fvMesh inletMesh
(
    io,
    pointField(infExtrude.points()),
    faceList(infExtrude.faces()),
    cellList(infExtrude.cells())
);

inletMesh.addFvPatches(infExtrude.patches(inletPatch));