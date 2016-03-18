/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  1.7.1                                 |
|   \\  /    A nd           | Web:      www.OpenFOAM.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    location    "0";
    object      Y.Liq;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 0 0 0 0 0 0];

internalField   uniform 0;

boundaryField
{
    movingWall
    {
        type            zeroGradient;
    }

    farFieldMoving
    {
        type            zeroGradient;
    }

    fixedWall
    {
        type            zeroGradient;
    }

    left
    {
        type            inletOutlet;
        value		uniform 0;
        inletValue	uniform 0;
    }

    farField
    {
        type            zeroGradient;
    }

    back
    {
        type            wedge;
    }

    front
    {
        type            wedge;
    }
    
    axis
    {
	type		empty;
    }
}

// ************************************************************************* //
