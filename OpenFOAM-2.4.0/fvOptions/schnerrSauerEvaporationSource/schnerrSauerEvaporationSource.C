/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "schnerrSauerEvaporationSource.H"
#include "fvMesh.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"
#include "basicThermo.H"
#include "coupledPolyPatch.H"
#include "surfaceInterpolate.H"
#include "fvm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(schnerrSauerEvaporationSource, 0);
    addToRunTimeSelectionTable
    (
        option,
        schnerrSauerEvaporationSource,
        dictionary
    );
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::schnerrSauerEvaporationSource::schnerrSauerEvaporationSource
(
    const word& name,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    option(name, modelType, dict, mesh),
    AlphaLiquidName_(coeffs_.subDict("fieldNames").lookup("liquidVolumeFraction")),
    YLiquidName_(coeffs_.subDict("fieldNames").lookup("liquidMassFraction")),
    pName_(coeffs_.subDict("fieldNames").lookup("pressure")),
    rhoName_(coeffs_.subDict("fieldNames").lookup("mixtureDensity")),
    rhoLiqName_(coeffs_.subDict("fieldNames").lookup("liquidDensity")),
    rhoGasName_(coeffs_.subDict("fieldNames").lookup("gasDensity")),
    pSat_(coeffs_.lookup("pSat")),
    n_(coeffs_.lookup("n")),
    dNuc_(coeffs_.lookup("dNuc")),
    Cc_(coeffs_.lookup("Cc")),
    Cv_(coeffs_.lookup("Cv"))
{
    fieldNames_.setSize(1, "YbarLiquid");
    applied_.setSize(1, false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::fv::schnerrSauerEvaporationSource::alphaNuc() const
{
    dimensionedScalar Vnuc = n_ * constant::mathematical::pi * pow3 (dNuc_) / 6;
    return Vnuc / (1.0 + Vnuc);
}

Foam::tmp<Foam::volScalarField>
Foam::fv::schnerrSauerEvaporationSource::rRb(const volScalarField& limitedAlphaLiq)
{
    return pow
    (
	((4 * constant::mathematical::pi * n_ ) / 3)
	*
	limitedAlphaLiq / (1.0 + alphaNuc() - limitedAlphaLiq),
	1.0 / 3.0
    );
}

bool Foam::fv::schnerrSauerEvaporationSource::alwaysApply() const
{
    return true;
}


void Foam::fv::schnerrSauerEvaporationSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label
)
{
    if (eqn.psi().name() != YLiquidName_)
    {
        return;
        #warning "raise exception here!"
    }
    
    const volScalarField& alphaLiq = mesh_.lookupObject<volScalarField>(AlphaLiquidName_);
    const volScalarField& p        = mesh_.lookupObject<volScalarField>(pName_);
    const volScalarField& rho      = mesh_.lookupObject<volScalarField>(rhoName_);
    const volScalarField& rhoLiq    = mesh_.lookupObject<volScalarField>(rhoLiqName_);
    const volScalarField& rhoGas    = mesh_.lookupObject<volScalarField>(rhoGasName_);

    const volScalarField limitedAlphaLiquid
    (
	"limitedAlphaLiquid",
	min(max(alphaLiq, 0.0), 1.0)
    );
    
    volScalarField limitedAlphaGas
    (
	"limitedAlphaGas",
	1.0 - limitedAlphaLiquid
    );
    
//    forAll(limitedAlphaGas, iCell)
//    {
//	if (limitedAlphaGas[iCell] < 1.0e-9)
//	{
//	    limitedAlphaGas[iCell] = 0.0;
//	}
//    }
    
    Info << "max/min of limitedAlphaLiquid: " << max(limitedAlphaLiquid).value() << "/" << min(limitedAlphaLiquid).value() << endl;
    Info << "max/min of limitedAlphaGas   : " << max(limitedAlphaGas).value() << "/" << min(limitedAlphaGas).value() << endl;

    dimensionedScalar p0 ("p0", p.dimensions(), 0.0);
    
    volScalarField pCoeff
    (
	"pCoeff",
	(rhoGas * rhoLiq / rho)
	* (3.0 * rRb(limitedAlphaLiquid))
	* sqrt
	(
	    (2. / 3.)
	    /
	    ((mag(p - pSat_) + 0.01*pSat_) * rhoLiq)
	)
    );
    
    volScalarField mLiqDotPlus
    (
	"mLiqDotPlus",
	Cc_ * limitedAlphaLiquid * limitedAlphaGas * pCoeff * max (p - pSat_, p0)
    );
    forAll(mLiqDotPlus.boundaryField(), iPatch)
    {
	forAll(mLiqDotPlus.boundaryField()[iPatch], iFace)
	{
	    mLiqDotPlus.boundaryField()[iPatch][iFace] = 0.0;
	}
    }
    Info << "Condensation source: " << max(mLiqDotPlus).value() << "/" << min(mLiqDotPlus).value() << endl;
    
    volScalarField mLiqDotMinus
    (
	"mLiqDotMinus",
	Cv_ * (rho / rhoLiq) * (limitedAlphaGas + alphaNuc()) * pCoeff * min(p - pSat_, p0)
    );
    forAll(mLiqDotMinus.boundaryField(), iPatch)
    {
	forAll(mLiqDotMinus.boundaryField()[iPatch], iFace)
	{
	    mLiqDotMinus.boundaryField()[iPatch][iFace] = 0.0;
	}
    }

    Info << "Evaporation source: " << max(mLiqDotMinus).value() << "/" << min(mLiqDotMinus).value() << endl;
    
//    eqn +=
//	(-mLiqDotPlus - fvm::Sp(mLiqDotMinus, eqn.psi()));
    eqn +=
    fvScalarMatrix
    (
	-fvm::Sp(mLiqDotMinus, eqn.psi())
	==
	mLiqDotPlus
    );
}


void Foam::fv::schnerrSauerEvaporationSource::writeData(Ostream& os) const
{
    os  << indent << name_ << endl;
    dict_.write(os);
}


bool Foam::fv::schnerrSauerEvaporationSource::read(const dictionary& dict)
{
    if (option::read(dict))
    {
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
