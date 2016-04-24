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
#include "fvc.H"

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
    AlphaLiquidName_(coeffs_.subDict("names").lookup("liquidVolumeFraction")),
    YLiquidName_(coeffs_.subDict("names").lookup("liquidMassFraction")),
    pName_(coeffs_.subDict("names").lookup("pressure")),
    rhoName_(coeffs_.subDict("names").lookup("mixtureDensity")),
    rhoLiqName_(coeffs_.subDict("names").lookup("liquidDensity")),
    rhoGasName_(coeffs_.subDict("names").lookup("gasDensity")),
    pSat_(coeffs_.lookup("pSat")),
    nMinus_(DataEntry<scalar>::New("nMinus", coeffs_)),
    nPlus_(DataEntry<scalar>::New("nPlus", coeffs_)),
    dNuc_(coeffs_.lookup("dNuc")),
    Cc_(coeffs_.lookup("Cc")),
    Cv_(coeffs_.lookup("Cv")),
    taur_(coeffs_.lookup("taur")),
    mLiqDotMinus_
    (
        IOobject
        (
            "mLiqDotMinus",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroMassSource", dimMass/dimLength/dimLength/dimLength/dimTime, 0)
    ),
    mLiqDotPlus_
    (
        IOobject
        (
            "mLiqDotPlus",
            mesh.time().timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("zeroMassSource", dimMass/dimLength/dimLength/dimLength/dimTime, 0)
    )
{
    fieldNames_.setSize(1, YLiquidName_);
    applied_.setSize(1, false);

    forAll(mLiqDotMinus_.boundaryField(), iPatch)
    {
	forAll(mLiqDotMinus_.boundaryField()[iPatch], iFace)
	{
	    mLiqDotMinus_.boundaryField()[iPatch][iFace] = 0.0;
	}
    }
    
    forAll(mLiqDotPlus_.boundaryField(), iPatch)
    {
	forAll(mLiqDotPlus_.boundaryField()[iPatch], iFace)
	{
	    mLiqDotPlus_.boundaryField()[iPatch][iFace] = 0.0;
	}
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dimensionedScalar
Foam::fv::schnerrSauerEvaporationSource::alphaNuc(const dimensionedScalar& nNuc) const
{
    dimensionedScalar Vnuc = nNuc * constant::mathematical::pi * pow3 (dNuc_) / 6;
    
    dimensionedScalar aNuc
    (
        "alphaNuc",
        dimless,
        max
        (
            SMALL,
            (Vnuc / (1.0 + Vnuc)).value()
        )
    );
    
    return aNuc;
}

Foam::tmp<Foam::volScalarField>
Foam::fv::schnerrSauerEvaporationSource::rRb(const volScalarField& limitedAlphaLiq, const dimensionedScalar& nNuc)
{
    return pow
    (
	((4 * constant::mathematical::pi * nNuc ) / 3)
	*
	limitedAlphaLiq / (1.0 + alphaNuc(nNuc) - limitedAlphaLiq),
	1.0 / 3.0
    );
}

bool Foam::fv::schnerrSauerEvaporationSource::alwaysApply() const
{
    return true;
}

void Foam::fv::schnerrSauerEvaporationSource::addSup
(
    const volScalarField& rhoMixture,
    fvMatrix<scalar>& eqn,
    const label fieldI
)
{
    const volScalarField& alphaLiq = mesh_.lookupObject<volScalarField>(AlphaLiquidName_);
    const volScalarField& p        = mesh_.lookupObject<volScalarField>(pName_);
    const volScalarField& rho      = mesh_.lookupObject<volScalarField>(rhoName_);
    const volScalarField& rhoLiq    = mesh_.lookupObject<volScalarField>(rhoLiqName_);
    const volScalarField& rhoGas    = mesh_.lookupObject<volScalarField>(rhoGasName_);

    const dimensionedScalar p0 ("p0", p.dimensions(), 0.0);

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
    
    mLiqDotMinus_.oldTime();
    mLiqDotPlus_.oldTime();
    
    dimensionedScalar nMinus
    (
        "nMinus",
        dimless / dimLength / dimLength / dimLength,
        nMinus_->value(mesh_.time().value())
    );
    
    dimensionedScalar nPlus
    (
        "nPlus",
        dimless / dimLength / dimLength / dimLength,
        nPlus_->value(mesh_.time().value())
    );
    
    volScalarField rhoRatio
    (
        rhoGas * rhoLiq / rho
    );
    
    volScalarField rhoByRbMinus
    (
	"rhoByRb",
	(rhoGas * rhoLiq / rho)
	* (3.0 * rRb(limitedAlphaLiquid, nMinus))
    );

    volScalarField rhoByRbPlus
    (
	"rhoByRb",
	(rhoGas * rhoLiq / rho)
	* (3.0 * rRb(limitedAlphaLiquid, nPlus))
    );

    volScalarField pCoeffMinus
    (
        "pCoeffMinus",
        sqrt(max(pSat_ - p, p0) / rhoLiq)
    );

    volScalarField pCoeffPlus
    (
        "pCoeffPlus",
        sqrt(max(p - pSat_, p0) / rhoLiq)
    );
    
    mLiqDotMinus_ = -Cv_ * (rho / rhoLiq)
        * (limitedAlphaGas + alphaNuc(nMinus)) * pCoeffMinus * rhoByRbMinus;
    mLiqDotPlus_  = Cc_ * limitedAlphaLiquid
        * limitedAlphaGas * pCoeffPlus * rhoByRbPlus;
    
    scalar rc = min(mesh_.time().deltaTValue() / taur_.value(), 1.0);
    
    mLiqDotMinus_.internalField() = mLiqDotMinus_.internalField()*rc + 
        mLiqDotMinus_.oldTime().internalField() * (1.0 - rc);
    mLiqDotPlus_.internalField() = mLiqDotPlus_.internalField()*rc + 
        mLiqDotPlus_.oldTime().internalField() * (1.0 - rc);
    
    Info<< "Evaporation source: "
        << gMin(mLiqDotMinus_.internalField()) << endl;
    Info<< "Condensation source: "
        << gMax(mLiqDotPlus_.internalField()) << endl;
    
    eqn.diag()   += (mLiqDotMinus_.internalField() * mLiqDotMinus_.mesh().V());
    eqn.source() += (-mLiqDotPlus_.internalField() * mLiqDotPlus_.mesh().V());
}



void Foam::fv::schnerrSauerEvaporationSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label
)
{
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
