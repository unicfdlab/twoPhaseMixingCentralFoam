{
    fvScalarMatrix rhoEqn
    (
        fvm::ddt(rho)
      + fvc::div(phi)
    );

    rhoEqn.solve();

    Info<< "rho min/max = " << min(rho).value() << "/"
        << max(rho).value() << endl;
}
