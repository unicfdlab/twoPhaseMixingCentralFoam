{
    word YiScheme ("div(phi," + YLiq.name() +")");
    
    surfaceScalarField rhoPhiYCorr
    (
        "rhoPhiYLiqCorr",
        fvc::flux
        (
            phi,
            YLiq,
            YiScheme
        )
    );
    
    surfaceScalarField mDiffC
    (
	"mDiffC",
	fvc::interpolate(turbulence->muEff()/Sc)
    );
    
    Info << "Creating matrix for fvOptions" << endl;
    
    #warning "change call to fvOptions to compressible"
    
    fvScalarMatrix SuSp
    (
	fvOptions(rho,YLiq)
    );
    
    mulesWithDiffusionImplicitLimiter
    (
	rho,
	YLiq,
	phi,
	allFacesLambda,
	rhoPhiYCorr,
	liqDiffusionFlux,
	mDiffC,
	SuSp
    );

    upwind<scalar> UDs(mesh, phi);

    fvScalarMatrix YLiqEqn
    (
        fvm::ddt(rho,YLiq)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phi, UDs).fvmDiv(phi, YLiq)
        +
        fvc::div(lambdaCoeffs * rhoPhiYCorr)
        ==
        SuSp
    );
    
    YLiqEqn.relax();
    
    YLiqEqn.solve(mesh.solver(YLiq.name()));
    
    surfaceScalarField centralSurfaceFlux
    (
        "centralSurfaceFlux",
        fvc::flux
        (
            phiPos,
            YLiq,
            YiScheme
        )
        +
        fvc::flux
        (
            phiNeg,
            YLiq,
            YiScheme
        )
    );

    upwind<scalar> UDsPos(mesh, phiPos);
    upwind<scalar> UDsNeg(mesh, phiNeg);
    
    surfaceScalarField centralSurfaceFluxBD
    (
        fv::gaussConvectionScheme<scalar>(mesh, phiPos, UDsPos).flux(phiPos, YLiq)
        +
        fv::gaussConvectionScheme<scalar>(mesh, phiNeg, UDsNeg).flux(phiNeg, YLiq)
    );
    
    //account for diffusion in enthalpy equation
    centralSurfaceFlux =
        centralSurfaceFluxBD + lambdaCoeffs*(centralSurfaceFlux - centralSurfaceFluxBD);

    liqDiffusionFlux = -
        (
            centralSurfaceFlux - (YLiqEqn.flux() + lambdaCoeffs * rhoPhiYCorr)
        );
    
    gasDiffusionFlux = -liqDiffusionFlux;

    Info<< "Mass fraction min/max = " << min(YLiq).value() << "/"
        << max(YLiq).value() << endl;
    
    YLiq = max(YLiq, 0.0);
    YLiq = min(YLiq, 1.0);
    YGas = 1.0 - YLiq;
    
}
