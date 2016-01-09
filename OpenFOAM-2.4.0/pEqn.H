
YbarLiq = YLiq * rhoGas / rhoLiq / (1.0 - YLiq + (rhoGas / rhoLiq) * YLiq);

Info << "max/min YbarLiq " << max(YbarLiq).value() << "/" << min(YbarLiq).value() << endl;

YbarLiq = min(max(YbarLiq,0.0),1.0);
YbarGas = 1.0 - YbarLiq;
//psi = YbarLiq * psiLiq + YbarGas * psiGas * (1.0 + YLiq * (rhoLiq0 / rhoLiq) / (1.0 - YLiq + YLiq * rhoGas / rhoLiq));
//alphaSqrRhoLiq0 = YbarLiq*YbarLiq*rhoLiq0;
psi = YbarLiq * psiLiq + YbarGas * psiGas  + (YbarGas * psiGas) * YLiq * (rhoLiq0 / rhoLiq) / (1.0 - YLiq + YLiq * rhoGas / rhoLiq);
alphaSqrRhoLiq0 = YbarLiq*YbarLiq*rhoLiq0;

const volScalarField p0 (p);

volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());
HbyA = UEqn.H() * rAU;

psi_pos = fvc::interpolate(psi, pos, "reconstruct(psi)");
psi_neg = fvc::interpolate(psi, neg, "reconstruct(psi)");

alphaSqrRhoLiq0_pos = fvc::interpolate(alphaSqrRhoLiq0, pos, "reconstruct(psi)");
alphaSqrRhoLiq0_neg = fvc::interpolate(alphaSqrRhoLiq0, neg, "reconstruct(psi)");

psiU_pos= fvc::interpolate(psi*HbyA, pos, "reconstruct(U)");
psiU_neg= fvc::interpolate(psi*HbyA, neg, "reconstruct(U)");

phiv_pos= (psiU_pos / psi_pos) & mesh.Sf();
phiv_neg= (psiU_neg / psi_neg) & mesh.Sf();

c = sqrt(thermo.Cp()/thermo.Cv() / psi);

//c = sqrt(1.0 / psi);

cSf_pos = fvc::interpolate(c, pos, "reconstruct(psi)")*mesh.magSf();
cSf_neg = fvc::interpolate(c, neg, "reconstruct(psi)")*mesh.magSf();

ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);

a_pos = ap/(ap - am);
aSf = am*a_pos;
a_neg = 1.0 - a_pos;

phiv_pos *= a_pos;
phiv_neg *= a_neg;
aphiv_pos = phiv_pos - aSf;
aphiv_neg = phiv_neg + aSf;
phid_pos = aphiv_pos * psi_pos;
phid_neg = aphiv_neg * psi_neg;

surfaceScalarField phiAlphaSqrRhoLiq0_pos
(
    "phiAlphaSqrRhoLiq0_pos",
    aphiv_pos*alphaSqrRhoLiq0_pos
);

surfaceScalarField phiAlphaSqrRhoLiq0_neg
(
    "phiAlphaSqrRhoLiq0_neg",
    aphiv_neg*alphaSqrRhoLiq0_neg
);

surfaceScalarField Dp_pos
(
    "Dp_pos",
    fvc::interpolate(rho*rAU, pos, "reconstruct(Dp)")
);

surfaceScalarField Dp_neg
(
    "Dp_neg",
    fvc::interpolate(rho*rAU, neg, "reconstruct(Dp)")
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn_pos
    (
	fvm::div(phid_pos,p) + fvc::div(phiAlphaSqrRhoLiq0_pos) - fvm::laplacian(Dp_pos*a_pos,p)
//	fvm::div(phid_pos,p) + fvc::div(phiAlphaSqrRhoLiq0_pos) - fvm::laplacian(rho*rAU,p)
    );
    
    fvScalarMatrix pEqn_neg
    (
	fvm::div(phid_neg,p) + fvc::div(phiAlphaSqrRhoLiq0_neg) - fvm::laplacian(Dp_neg*a_neg,p)
	//fvm::div(phid_neg,p) + fvc::div(phiAlphaSqrRhoLiq0_neg)
    );
    
    fvScalarMatrix pEqn
    (
	fvm::ddt(psi,p)
	+
	fvc::ddt(alphaSqrRhoLiq0)
	+
	pEqn_pos
	+
	pEqn_neg
    );
    
    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
    if (pimple.finalNonOrthogonalIter())
    {
	phiPos = pEqn_pos.flux() + phiAlphaSqrRhoLiq0_pos;
	phiNeg = pEqn_neg.flux() + phiAlphaSqrRhoLiq0_neg;
	phi = phiPos + phiNeg;
    }
}

p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
p_neg = fvc::interpolate(p, neg, "reconstruct(p)");

gradp = fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf());

#include "limitPressureGradient.H"
U = HbyA - rAU * gradp;
U.correctBoundaryConditions();

K = 0.5*magSqr(U);
    
// Update EOS
rhoLiq += psiLiq * (p - p0);
rhoGas += psiGas * (p - p0);

Info << "max(mag(U)): " << max(mag(U)).value() << endl;

solve
(
    fvm::ddt(rho) + fvc::div(phi)
);

volScalarField rhoT = 1.0 / (YLiq / rhoLiq + YGas / rhoGas);

#include "compressibleContinuityErrs.H"

rho = rhoT;


