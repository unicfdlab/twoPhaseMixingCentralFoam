{
    word YiPosScheme ("div(phiPos," + YLiq.name() +")");
    word YiNegScheme ("div(phiNeg," + YLiq.name() +")");

    surfaceScalarField rhoPhiYCorr
    (
	"rhoPhiYLiqCorr",
	
	fvc::flux
	(
	    phiPos,
	    YLiq,
	    YiPosScheme
	)
	+
	fvc::flux
	(
	    phiNeg,
	    YLiq,
	    YiNegScheme
	)
    );
    
    surfaceScalarField mDiffC
    (
	"mDiffC",
	fvc::interpolate(turbulence->muEff()/Sc)
    );
    
    fvScalarMatrix SuSp
    (
	fvOptions(rho, YLiq)
    );
    
    mulesWithDiffusionImplicitLimiter
    (
	rho,
	YLiq,
	phiPos,
	phiNeg,
	allFacesLambda,
	rhoPhiYCorr,
	liqDiffusionFlux,
	mDiffC,
	SuSp
    );
    
    fvScalarMatrix YLiqDiffusion
    (
	- fvm::laplacian(mDiffC, YLiq)
    );

    upwind<scalar> UDsPos(mesh, phiPos);
    upwind<scalar> UDsNeg(mesh, phiNeg);

    fvScalarMatrix YLiqConvection
    (
	fvm::ddt(rho,YLiq)
	+
	fv::gaussConvectionScheme<scalar>(mesh, phiPos, UDsPos).fvmDiv(phiPos, YLiq)
	+
	fv::gaussConvectionScheme<scalar>(mesh, phiNeg, UDsNeg).fvmDiv(phiNeg, YLiq)
	==
	SuSp
    );
    
    fvScalarMatrix YLiqEqn
    (
	YLiqConvection + YLiqDiffusion +  fvc::div(lambdaCoeffs * rhoPhiYCorr)
    );
    
    YLiqEqn.relax();
    
    YLiqEqn.solve(mesh.solver(YLiq.name()));
    
    Info<< "Mass fraction min/max = " << min(YLiq).value() << "/"
        << max(YLiq).value() << endl;
    
    YLiq = max(YLiq, 0.0);
    YLiq = min(YLiq, 1.0);
    YGas = 1.0 - YLiq;
        
    liqDiffusionFlux = YLiqDiffusion.flux();
    gasDiffusionFlux = -liqDiffusionFlux;
}
