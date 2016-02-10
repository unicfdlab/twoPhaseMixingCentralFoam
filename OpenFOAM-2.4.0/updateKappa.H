{
    aphiv_pos = phiPos / (p_pos * psi_pos + rhoHat_pos);
    aphiv_neg = phiNeg / (p_neg * psi_neg + rhoHat_neg);
    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

    surfaceScalarField amaxSfbyDelta
    (
        mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
    );
    
    surfaceScalarField Maf
    (
        mag(phi) / (psi_pos*p_pos*a_pos + psi_neg*p_neg*a_neg + rhoHat_pos*a_pos + rhoHat_neg*a_neg)
        / (cSf_pos*a_pos + cSf_neg*a_neg)
    );

    Info << "max/min Maf: " << max(Maf).value() << "/" << min(Maf).value() << endl;

    const scalar MaTreshold
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("MaThreshold")
        )
    );

    const scalar Ampl
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("Amplitude")
        )
    );
    
    const scalar Slope
    (
        readScalar
        (
            mesh.solutionDict().subDict("PIMPLE").lookup("Slope")
        )
    );

    surfaceScalarField MaFunction
    (
        "MaFunction",
        Maf*Ampl*(tanh(Slope*(Maf - MaTreshold)) + 1)
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

