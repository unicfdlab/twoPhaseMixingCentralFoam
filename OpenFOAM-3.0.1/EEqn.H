{
    volScalarField& he = thermo.he();
    volScalarField alphahEff (turbulence->alphaEff());
    volScalarField Cp = thermo.Cp();

    upwind<scalar> UDsPos(mesh, phiPos);
    upwind<scalar> UDsNeg(mesh, phiNeg);
    upwind<scalar> UDs(mesh, phi);

    fvScalarMatrix hConvection
    (
        fvm::ddt(rho,he)
        //+
        //fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, he)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phiPos, UDsPos).fvmDiv(phiPos, he)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phiNeg, UDsNeg).fvmDiv(phiNeg, he)
    );

    /*
     *
     * Explicit flux correction, bounded with MULES
     *
     */
    
    surfaceScalarField hFluxCorr
    (
        lambdaCoeffs
        *
        (
	    fvc::flux
	    (
		phiPos,
		he,
		"div(phiPos," + he.name() + ")"
	    )
	    +
	    fvc::flux
	    (
		phiNeg,
		he,
		"div(phiNeg," + he.name() + ")"
	    )
	    -
	    hConvection.flux()
        )
    );
    
    surfaceScalarField hCmptFlux
    (
        "hCmptFlux",
        fvc::flux
        (
            liqDiffusionFlux,
            thermo.thermoLiq().he(),
            "div(rhoi*Uri,hi)"
        )
        +
        fvc::flux
        (
            gasDiffusionFlux,
            thermo.thermoGas().he(),
            "div(rhoi*Uri,hi)"
        )
        +
        fvc::flux
        (
            fvc::snGrad(YLiq) * mesh.magSf(),
            thermo.thermoLiq().he() * alphahEff,
            "div(alphaEff*hi,grad(Yi))"
        )
        +
        fvc::flux
        (
            fvc::snGrad(YGas) * mesh.magSf(),
            thermo.thermoGas().he() * alphahEff,
            "div(alphaEff*hi,grad(Yi))"
        )
    );
    
    fvScalarMatrix EEqn
    (
        hConvection + fvc::div(hFluxCorr) + fvc::div(hCmptFlux) 
        +
        EkChange
        ==
        fvm::laplacian(alphahEff, he)
        +
        dpdt
//        +
//        rho*(g & U)
//        +
//        fvOptions(rho,he)
    );
    
    fvOptions.constrain(EEqn);
    
    EEqn.solve();
    
    fvOptions.correct(he);
    
    thermo.correct();
    
    //surfaceScalarField eFlux = EEqn.flux() + hCmptFlux - hConvection.flux() - hFluxCorr;
    //
    //Info << "Max/min of eFlux: " << max(eFlux).value() << "/" << min(eFlux).value() << endl;
    //he.write();
    //alphahEff.rename("alphahEff");
    //alphahEff.write();
    
    //T.correctBoundaryConditions();
    
}
