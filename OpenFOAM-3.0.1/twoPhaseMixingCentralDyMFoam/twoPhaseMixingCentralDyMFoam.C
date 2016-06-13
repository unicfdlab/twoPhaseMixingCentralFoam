/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    twoPhaseMixingCentralDyMFoam

Description
    Transient Eulerian two-phase solver with dynamic meshes. Liquid and gas are
    considered as compressible fluids. Mass transfer at the interface
    is accounted at the diffusion approximation.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "pimpleControl.H"
#include "kappaFunction.H"
#include "customMULES.H"
#include "fvIOoptionList.H"
#include "cellQuality.H"
#include "fvcSmooth.H"
#include "compressibleTwoPhaseMixtureThermo.H"
#include "dynamicFvMesh.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    pimpleControl pimple(mesh);
    
    #include "createFields.H"
    #include "createMRF.H"
    #include "createFvOptions.H"
    bool checkMeshCourantNo =
                readBool(pimple.dict().lookup("checkMeshCourantNo"));

    #include "createTimeControls.H"
    #include "readTimeControls.H"

    dimensionedScalar v_zero ("v_zero", dimVolume/dimTime, 0.0);
    #include "createSurfaceFields.H"
    #include "markBadQualityCells.H"
    #include "readCourantType.H"
    #include "centralCompressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    surfaceScalarField meshPhi
    (
        "volMeshPhi",
        phiv_pos * 0.0
    );

    surfaceScalarField meshPhiPos
    (
        "meshPhiPos",
        meshPhi * rho_pos * 0.0
    );

    surfaceScalarField meshPhiNeg
    (
        "meshPhiPos",
        meshPhi * rho_neg * 0.0
    );

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
	#include "acousticCourantNo.H"
        #include "centralCompressibleCourantNo.H"
        #include "setDeltaT.H"

        runTime++;
        
        thermo.he().oldTime();
        rho.oldTime();
        p.oldTime();
        psi.oldTime();
        YLiq.oldTime();
        YGas.oldTime();
        YbarLiq.oldTime();
        YbarGas.oldTime();
        rhoHat.oldTime();
        K.oldTime();

        Info<< "Time = " << runTime.timeName() << nl << endl;
        // --- Move mesh and update fluxes
        {
        
            // Do any mesh changes
            mesh.update();
          
             if (mesh.changing())
             {
                 meshPhi = fvc::meshPhi(rho,U)();
                 
                 meshPhiPos = a_pos*rho_pos*meshPhi;
                 meshPhiNeg = a_neg*rho_neg*meshPhi;
                 //make fluxes relative
                 phiPos -= (meshPhiPos + (1.0 - kappa)*meshPhiNeg);
                 phiNeg -= (meshPhiNeg*kappa);
                 phi = phiPos + phiNeg;

                 if (checkMeshCourantNo)
                 {
                     #include "customMeshCourantNo.H"
                 }

                 #include "markBadQualityCells.H"
             }
        }

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "MixtureRhoEqn.H"
        
            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        
            scalarField allFacesLambda(mesh.nFaces(), 1.0);
            slicedSurfaceScalarField lambdaCoeffs
            (
                IOobject
                (
                    "lambdaCoeffs",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh,
                dimless,
                allFacesLambda,
                false   // Use slices for the couples
            );
            
            #include "YLiqEqn.H"
            #include "UEqn.H"
            #include "EEqn.H"
            
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqnDyM.H"
            }
            
            K = 0.5*magSqr(U);
            dpdt = (p - p.oldTime()) / runTime.deltaT();
            EkChange = fvc::ddt(rho,K) + fvc::div(phi,K)
                - fvc::div( ((-turbulence->devRhoReff()) & U) );
            //EkChange = fvc::ddt(rho,K) + fvc::div(phiPos,K) + fvc::div(phiNeg,K)
            //  - fvc::div( ((-turbulence->devRhoReff()) & U) );
            
            //make fluxes absolute
            phiPos += meshPhiPos;
            phiNeg += meshPhiNeg;
            phi = phiPos + phiNeg;
            #include "updateKappa.H"
        }
        
        if(runTime.write())
        {
            c.write();
            YbarLiq.write();
        }
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
