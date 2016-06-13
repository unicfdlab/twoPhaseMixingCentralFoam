{
    aphiv_pos = phiPos / rho_pos;
    aphiv_neg = phiNeg / rho_neg;
    dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(mesh.magSf()).value());
    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg)) + amaxSmall;
    
    kappaFuncPtr().update();
    kappa = kappaFuncPtr().kappa()();
//    surfaceScalarField MaFunction
//    (
//        "MaFunction",
//        //Maf*Ampl*(tanh(Slope*(Maf - MaTreshold)) + 1)
//        //Maf*0.0 + 1.0
//        //Maf
//        max(Maf,min(FaceAcCourant,1.0))
//    );
    
//    Info << "max/min MaFunction: " << max(MaFunction).value() << "/" << min(MaFunction).value() << endl;

    Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
    phiPos = phiPos + (1.0 - kappa) * phiNeg;
    phiNeg = kappa * phiNeg;
    
    //phiNeg += (1.0 - kappa) * phiPos;
    //phiPos *= kappa;
}

