
if (pimple.corrPISO() > 1)
{
    thermo.correct();
}

volScalarField rhoT = thermo.rho();
rhoHat = rhoT - psi*p;

volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());

HbyA = UEqn.H() * rAU;
{
    rho_pos     = fvc::interpolate(rho, pos, "reconstruct(rho)");
    rho_neg     = fvc::interpolate(rho, neg, "reconstruct(rho)");

    psi_pos     = fvc::interpolate(psi, pos, "reconstruct(psi)");
    psi_neg     = fvc::interpolate(psi, neg, "reconstruct(psi)");

    rhoHat_pos  = fvc::interpolate(rhoHat, pos, "reconstruct(rhoHat)");
    rhoHat_neg  = fvc::interpolate(rhoHat, neg, "reconstruct(rhoHat)");

    rhoU_pos    = fvc::interpolate(rho*HbyA, pos, "reconstruct(U)");
    rhoU_neg    = fvc::interpolate(rho*HbyA, neg, "reconstruct(U)");

    phiv_pos    = (rhoU_pos / rho_pos) & mesh.Sf();
    phiv_neg    = (rhoU_neg / rho_neg) & mesh.Sf();

    c           = sqrt(thermo.Cp()/thermo.Cv() / psi);
    cSf_pos     = fvc::interpolate(c, pos, "reconstruct(c)")*mesh.magSf();
    cSf_neg     = fvc::interpolate(c, neg, "reconstruct(c)")*mesh.magSf();

    ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
    am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);
    
    YbarLiq_pos = fvc::interpolate(YbarLiq, pos, "reconstruct(YbarLiq)");
    YbarLiq_neg = fvc::interpolate(YbarLiq, neg, "reconstruct(YbarLiq)");

    a_pos       = ap/(ap - am);
    aSf         = am*a_pos;
    a_neg       = 1.0 - a_pos;

    aphiv_pos   = a_pos*phiv_pos - aSf;
    aphiv_neg   = a_neg*phiv_neg + aSf;

    phid_pos = aphiv_pos * psi_pos;
    phid_neg = aphiv_neg * psi_neg;

    phiRhoHat_pos = aphiv_pos * rhoHat_pos;
    phiRhoHat_neg = aphiv_neg * rhoHat_neg;

    //update anti-diffusion fluxes in the interface regions
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
            phiPos = pEqn_pos.flux() + phiRhoHat_pos;
            phiNeg = pEqn_neg.flux() + phiRhoHat_neg;
            phi = phiPos + phiNeg;
        }
    }

    p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
    p_neg = fvc::interpolate(p, neg, "reconstruct(p)");
    rho_pos = p_pos*psi_pos + rhoHat_pos;
    rho_neg = p_neg*psi_neg + rhoHat_neg;

    gradp = fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf());

    #include "limitPressureGradient.H"
}

U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

K = 0.5*magSqr(U);

Info << "max(mag(U)): " << max(mag(U)).value() << endl;
Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

#include "MixtureRhoEqn.H"

p = max(p, pMin);
p.correctBoundaryConditions();
rhoT = rhoHat + psi*p; //1.0 / (YLiq / rhoLiq + YGas / rhoGas);

#include "compressibleContinuityErrs.H"

rho = rhoT;


