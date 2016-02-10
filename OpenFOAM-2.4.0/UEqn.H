    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phiPos, U)
      + fvm::div(phiNeg, U)
      + fvm::div(liqDiffusionFlux, U)
      + fvm::div(gasDiffusionFlux, U)
      + turbulence->divDevRhoReff(U)
     ==
        rho.dimensionedInternalField()*g
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -gradp);

        K = 0.5*magSqr(U);
    }
