{
    volScalarField& he = thermo.he();
    volScalarField alphaEff (turbulence->alphaEff());
    volScalarField Cp = thermo.Cp();

    upwind<scalar> UDsPos(mesh, phiPos);
    upwind<scalar> UDsNeg(mesh, phiNeg);

    fvScalarMatrix hConvection
    (
	fvm::ddt(rho,he)
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
    );

    fvScalarMatrix EEqn
    (
	hConvection + fvc::div(hFluxCorr) + fvc::div(hCmptFlux)
	+
	EkChange
	==
	fvm::laplacian(alphaEff, he)
	-
	fvc::laplacian(alphaEff*T, thermo.Cp())
	+
	dpdt
	+
	rho*(g & U)
    );

    EEqn.solve();
    
    thermo.correct();
    
    T.correctBoundaryConditions();
}
