    surfaceScalarField rho_pos
    (
        "rho_pos",
        fvc::interpolate(thermo.rho(), pos, "reconstruct(rho)")
    );
    surfaceScalarField rho_neg
    (
        "rho_neg",
        fvc::interpolate(thermo.rho(), neg, "reconstruct(rho)")
    );

    surfaceScalarField psi_pos
    (
        "psi_pos",
        fvc::interpolate(psi, pos, "reconstruct(psi)")
    );
    surfaceScalarField psi_neg
    (
        "psi_neg",
        fvc::interpolate(psi, neg, "reconstruct(psi)")
    );

    surfaceVectorField rhoU_pos
    (
        "rhoU_pos",
        fvc::interpolate(rho*U, pos, "reconstruct(U)")
    );
    surfaceVectorField rhoU_neg
    (
        "rhoU_neg",
        fvc::interpolate(rho*U, neg, "reconstruct(U)")
    );

    surfaceScalarField p_pos
    (
	"p_pos",
	fvc::interpolate(p, pos, "reconstruct(p)")
    );
    surfaceScalarField p_neg
    (
	"p_neg",
	fvc::interpolate(p, neg, "reconstruct(p)")
    );
    
    surfaceScalarField phiv_pos
    (
	"phiv_pos",
	(rhoU_pos / rho_pos)  & mesh.Sf()
    );
    surfaceScalarField phiv_neg
    (
	"phiv_neg",
	(rhoU_neg / rho_neg) & mesh.Sf()
    );

    volScalarField c
    (
	"cSound",
	sqrt(thermo.Cp()/thermo.Cv() / psi)
    );
    surfaceScalarField cSf_pos
    (
	"cSf_pos",
	fvc::interpolate(c, pos, "reconstruct(psi)")*mesh.magSf()
    );
    surfaceScalarField cSf_neg
    (
	"cSf_neg",
	fvc::interpolate(c, neg, "reconstruct(psi)")*mesh.magSf()
    );
    
    surfaceScalarField ap
    (
	"ap",
	max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
    );
    surfaceScalarField am
    (
	"am",
	min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
    );

    surfaceScalarField a_pos
    (
	"a_pos",
	ap/(ap - am)
    );
    surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));
    surfaceScalarField aSf
    (
	"aSf",
	am*a_pos*0.0
    );
    surfaceScalarField a_neg
    (
	"a_neg",
	1.0 - a_pos
    );
    
    phiv_pos *= a_pos;
    phiv_neg *= a_neg;
    
    surfaceScalarField aphiv_pos("aphiv_pos", phiv_pos - aSf);
    surfaceScalarField aphiv_neg("aphiv_neg", phiv_neg + aSf);
    
    surfaceScalarField phid_pos
    (
	"phid_pos",
	aphiv_pos * psi_pos
    );
    surfaceScalarField phid_neg
    (
	"phid_neg",
	aphiv_neg * psi_neg
    );

    surfaceScalarField rhoHat_pos
    (
        "rhoHat_pos",
        fvc::interpolate(rhoHat, pos, "reconstruct(rhoHat)")
    );
    
    surfaceScalarField rhoHat_neg
    (
        "rhoHat_neg",
        fvc::interpolate(rhoHat, neg, "reconstruct(rhoHat)")
    );
    
    surfaceScalarField liqDiffusionFlux
    (
        "liqDiffusionFlux",
        fvc::interpolate(turbulence->muEff()) * ( (fvc::snGrad(YLiq)) * mesh.magSf()) * 0.0
    );
    
    surfaceScalarField gasDiffusionFlux
    (
        "gasDiffusionFlux",
        liqDiffusionFlux
    );

    surfaceScalarField phiPos
    (
	"phiPos",
	phi
    );
    surfaceScalarField phiNeg
    (
	"phiNeg",
	phi * 0.0
    );

    surfaceScalarField phiRhoHat_pos
    (
        "phiRhoHat_pos",
        phiPos * 0.0
    );
    surfaceScalarField phiRhoHat_neg
    (
        "phiRhoHat_neg",
        phiRhoHat_pos
    );

    surfaceScalarField YbarLiq_pos
    (
        "YbarLiq_pos",
        fvc::interpolate(YbarLiq, pos, "reconstruct(YbarLiq)")
    );
    surfaceScalarField YbarLiq_neg
    (
        "YbarLiq_neg",
        fvc::interpolate(YbarLiq, neg, "reconstruct(YbarLiq)")
    );
    
    surfaceScalarField Dp_pos
    (
        "Dp_pos",
        a_pos * dimensionedScalar ("zeroA", dimTime, 0.0)
    );
    
    surfaceScalarField Dp_neg
    (
        "Dp_neg",
        Dp_pos
    );
    
    /*
     *
     * Flux scheme switching surface field
     *
     */
    
    autoPtr<fv::kappaFunction> kappaFuncPtr
    (
        fv::kappaFunction::New
        (
            "kappaFunction",
            mesh.solutionDict().subDict("PIMPLE").subDict("kappaFunction"),
            mesh
        )
    );
    
    surfaceScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    );
    
    forAll(kappa, iFace)
    {
        kappa[iFace] = 0.0;
    }
    
    forAll(kappa.boundaryField(), iPatch)
    {
        kappa.boundaryField()[iPatch] = 0.0;
    }

