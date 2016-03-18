    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phiPos, U)
      + fvm::div(phiNeg, U)
      + fvm::div(liqDiffusionFlux, U)
      + fvm::div(gasDiffusionFlux, U)
      + turbulence->divDevRhoReff(U)
      + MRF.DDt(rho,U)
     ==
        rho.dimensionedInternalField()*g
    );

    fvOptions.constrain(UEqn);
    
    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn == -gradp
        );
        
        fvOptions.correct(U);
        
        K = 0.5*magSqr(U);
    }
    
    UEqn.relax();


