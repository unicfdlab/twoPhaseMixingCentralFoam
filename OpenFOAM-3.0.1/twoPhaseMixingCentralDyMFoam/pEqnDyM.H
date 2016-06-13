
if (pimple.corrPISO() > 1)
{
    thermo.correct();
}

volScalarField rhoT = thermo.rho();
rhoHat = rhoT - psi*p;
volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());

HbyA = UEqn.H() * rAU;

/*
 *
 * Standard FVM p-U coupling procedure
 *
 */
if(false)
{
    surfaceScalarField psif
    (
        "psif",
        fvc::interpolate(psi)
    );

    surfaceScalarField phiv
    (
        "phiv",
        (fvc::interpolate(HbyA) & mesh.Sf()) - meshPhi
    );

    surfaceScalarField phid
    (
        "phid",
        phiv * psif
    );

    surfaceScalarField phiHat
    (
        "phiHat",
        phiv * fvc::interpolate(rhoHat)
    );

    surfaceScalarField Dp
    (
        "Dp",
        fvc::interpolate(rho*rAU)
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
            +
            fvc::ddt(rhoHat)
            +
            fvc::div(phiHat)
            +
            fvm::div(phid,p)
            -
            fvm::laplacian(Dp,p)
            ==
            fvOptions(psi,p,rho.name())
        );
    
        fvOptions.constrain(pEqn);
    
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
        if (pimple.finalNonOrthogonalIter())
        {
	    phi = phiHat + pEqn.flux();
        }
    }
    
    gradp = fvc::grad(p);
}


/*
 *
 * Hybrid FVM p-U coupling procedure
 *
 */
if(true)
{
    rho_pos       = fvc::interpolate(rho, pos, "reconstruct(rho)");
    rho_neg       = fvc::interpolate(rho, neg, "reconstruct(rho)");
    
    psi_pos       = fvc::interpolate(psi, pos, "reconstruct(psi)");
    psi_neg       = fvc::interpolate(psi, neg, "reconstruct(psi)");
    
    rhoHat_pos    = fvc::interpolate(rhoHat, pos, "reconstruct(rhoHat)");
    rhoHat_neg    = fvc::interpolate(rhoHat, neg, "reconstruct(rhoHat)");
    
    rhoU_pos      = fvc::interpolate(rho*HbyA, pos, "reconstruct(U)");
    rhoU_neg      = fvc::interpolate(rho*HbyA, neg, "reconstruct(U)");

    phiv_pos      = (rhoU_pos / rho_pos) & mesh.Sf();
    phiv_neg      = (rhoU_neg / rho_neg) & mesh.Sf();
    
    c             = sqrt(thermo.Cp()/thermo.Cv() / psi);
    cSf_pos       = fvc::interpolate(c, pos, "reconstruct(c)")*mesh.magSf();
    cSf_neg       = fvc::interpolate(c, neg, "reconstruct(c)")*mesh.magSf();
    
    phiv_pos -= meshPhi;
    phiv_neg -= meshPhi;
    
    MRF.makeRelative(phiv_pos);
    MRF.makeRelative(phiv_neg);
    
    ap            = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
    am            = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);
    
    a_pos         = ap/(ap - am);
    aSf           = am*a_pos;
    a_neg         = 1.0 - a_pos;
    
//    aphiv_pos     = a_pos * phiv_pos + aSf;
//    aphiv_neg     = a_neg * phiv_neg - aSf;
    
    aphiv_pos     = a_pos * phiv_pos - aSf;
    aphiv_neg     = a_neg * phiv_neg + aSf;
    
    phid_pos      = aphiv_pos * psi_pos;
    phid_neg      = aphiv_neg * psi_neg;
    
    phiRhoHat_pos = aphiv_pos * rhoHat_pos;
    phiRhoHat_neg = aphiv_neg * rhoHat_neg;
    
    YbarLiq_pos   = fvc::interpolate(YbarLiq, pos, "reconstruct(rhoHat)");
    YbarLiq_neg   = fvc::interpolate(YbarLiq, neg, "reconstruct(rhoHat)");
    
    surfaceScalarField meshPhid_pos
    (
        "meshPhid",
        a_pos*psi_pos*meshPhi
    );
    surfaceScalarField meshPhid_neg
    (
        "meshPhid",
        a_neg*psi_neg*meshPhi
    );
    surfaceScalarField meshPhiRhoHat_pos
    (
        "meshPhiRhoHat_pos",
        a_pos*rhoHat_pos*meshPhi
    );
    surfaceScalarField meshPhiRhoHat_neg
    (
        "meshPhiRhoHat_neg",
        a_neg*rhoHat_neg*meshPhi
    );

    //remove diffusion term for faces where mass fraction is changing (Cid == 0 at fvSolution)
    {
        phiRhoHat_pos -=
            aSf * (YbarLiq_neg - YbarLiq_pos) * rhoHat_pos * Cid;
    
        phiRhoHat_neg -=
            aSf * (YbarLiq_neg - YbarLiq_pos) * rhoHat_neg * Cid;
    }

    Dp_pos = a_pos*fvc::interpolate(rho*rAU, pos, "reconstruct(Dp)");

    Dp_neg = a_neg*fvc::interpolate(rho*rAU, neg, "reconstruct(Dp)");

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn_pos
        (
            fvm::div(phid_pos,p) + fvc::div(phiRhoHat_pos) - fvm::laplacian(Dp_pos,p)
        );
    
        fvScalarMatrix pEqn_neg
        (
            fvm::div(phid_neg,p) + fvc::div(phiRhoHat_neg) - fvm::laplacian(Dp_neg,p)
        );
    
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi,p)
            +
            fvc::ddt(rhoHat)
            +
            pEqn_pos
            +
            pEqn_neg
            ==
            fvOptions(psi,p,rho.name())
        );
    
        fvOptions.constrain(pEqn);
    
        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
        if (pimple.finalNonOrthogonalIter())
        {
            phiPos = pEqn_pos.flux() + phiRhoHat_pos;//Relative
            phiNeg = pEqn_neg.flux() + phiRhoHat_neg;//Relative
            phi = phiPos + phiNeg;                   //Relative
        }
    }
    
    p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
    p_neg = fvc::interpolate(p, neg, "reconstruct(p)");
    
    gradp = fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf());

    #include "limitPressureGradient.H"
}

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

Info << "max(mag(U)): " << max(mag(U)).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

//rhoT = rhoHat + psi*p; //1.0 / (YLiq / rhoLiq + YGas / rhoGas);
thermo.correctRealDensities();
rhoT = 1.0 / (YLiq / thermo.thermoLiq().rho() + YGas / thermo.thermoGas().rho());
#include "MixtureRhoEqn.H"

#include "compressibleContinuityErrs.H"
rho = rhoT;
Info << "max/min rhoT: " << max(rho).value() << "/" << min(rhoT).value() << endl;

rho_pos = fvc::interpolate(rho, pos, "reconstruct(rho)");
rho_neg = fvc::interpolate(rho, neg, "reconstruct(rho)");

meshPhiPos = a_pos*rho_pos*meshPhi;
meshPhiNeg = a_neg*rho_neg*meshPhi;

//
//END-OF-FILE
//

