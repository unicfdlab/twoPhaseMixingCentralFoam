{
    aphiv_pos = phiPos / rho_pos;
    aphiv_neg = phiNeg / rho_neg;
    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

    surfaceScalarField amaxSfbyDelta
    (
        mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
    );
    
    surfaceScalarField Maf
    (
        max
        (
            mag(phi) / (rho_pos*a_pos + rho_neg*a_neg)
            / (cSf_pos*a_pos + cSf_neg*a_neg),
            scalar(0)
        )
    );

    Info << "max/min Maf: " << max(Maf).value() << "/" << min(Maf).value() << endl;

//    const scalar MaTreshold
//    (
//        readScalar
//        (
//            mesh.solutionDict().subDict("PIMPLE").lookup("MaThreshold")
//        )
//    );
//
//    const scalar Ampl
//    (
//        readScalar
//        (
//            mesh.solutionDict().subDict("PIMPLE").lookup("Amplitude")
//        )
//    );
//    
//    const scalar Slope
//    (
//        readScalar
//        (
//            mesh.solutionDict().subDict("PIMPLE").lookup("Slope")
//        )
//    );

    surfaceScalarField MaFunction
    (
        "MaFunction",
        //Maf*Ampl*(tanh(Slope*(Maf - MaTreshold)) + 1)
        Maf
    );
    
    Info << "max/min MaFunction: " << max(MaFunction).value() << "/" << min(MaFunction).value() << endl;
    
    surfaceScalarField FaceAcCourant
    (
        "FaceAcCourant",
        (amaxSfbyDelta/mesh.magSf() * runTime.deltaT())
    );

    Info << "max/min FaceAcCourant: " << max(FaceAcCourant).value() << "/" << min(FaceAcCourant).value() << endl;
    
    kappa = 
        min
        (
            MaFunction / FaceAcCourant,
            scalar(1.0)
        );

    forAll(kappa.boundaryField(), iPatch)
    {
        fvsPatchField<scalar>& kappapf = kappa.boundaryField()[iPatch];
        if (isA<coupledFvsPatchField<scalar> > (kappapf))
        {
            forAll (kappapf, iFace)
            {
                kappapf[iFace] = 0.0;
            }
        }
    }

    Info << "max/min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
    phiPos = phiPos + (1.0 - kappa) * phiNeg;
    phiNeg = kappa * phiNeg;
}

