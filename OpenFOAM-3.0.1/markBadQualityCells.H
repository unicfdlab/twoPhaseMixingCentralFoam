{
    Info << "Marking cells with bad quality to limit pressure gradient" << endl;
    
    badQualityCells.resize(0);
    
    const scalar nonOrthogonalityThreshold
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("nonOrthogonalityThreshold")
        )
    );
    
    const scalar skewnessThreshold
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("skewnessThreshold")
        )
    );

//    volScalarField cellMarker
//    (
//        "cellMarker",
//        YLiq * 0.0
//    );

    cellQuality cq (mesh);
    
    scalarField cellNonOrtho (cq.nonOrthogonality());
    
    scalarField cellSkewness (cq.skewness());
    
    forAll(cellNonOrtho, iCell)
    {
        if ( (cellNonOrtho[iCell] >= nonOrthogonalityThreshold) ||
            (cellSkewness[iCell] >= skewnessThreshold) )
        {
            badQualityCells.append(iCell);
            
//            cellMarker[iCell] = 1.0;
        }
    }
    
//    if (runTime.outputTime())
//    {
//        cellMarker.write();
//    }
}
