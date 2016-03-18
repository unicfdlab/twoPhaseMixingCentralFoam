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

Class
    Foam::fv::schnerrSauerEvaporationSource

Description
    Schnerr-Sauer cavitation model, based on reducted Rayleigh-Plesset
    equation




SourceFiles
    schnerrSauerEvaporationSource.C

\*---------------------------------------------------------------------------*/

#ifndef schnerrSauerEvaporationSource_H
#define schnerrSauerEvaporationSource_H

#include "fvOption.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"
#include "DataEntry.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
              Class schnerrSauerEvaporationSource Declaration
\*---------------------------------------------------------------------------*/

class schnerrSauerEvaporationSource
:
    public option
{

protected:

    // Protected data
	
        //- Name of liquid volume fraction field
        word AlphaLiquidName_;
        
        //- Name of liquid mass fraction field
        word YLiquidName_;
	
        //- Name of pressure field
        word pName_;
        
        //- Name of mixture density field
        word rhoName_;
        
        //-
        word rhoLiqName_;
        
        //-
        word rhoGasName_;
        
        //- Saturation pressure
        dimensionedScalar pSat_;
        
        //- Bubble number density
        autoPtr<DataEntry<scalar> > nMinus_;

        //- Bubble number density
        autoPtr<DataEntry<scalar> > nPlus_;

        //- Nucleation site diameter
        dimensionedScalar dNuc_;
        
        //- Condensation constant
        dimensionedScalar Cc_;
        
        //- Vaporization constant
        dimensionedScalar Cv_;
        
        //-
        dimensionedScalar taur_;
        
        //-
        volScalarField mLiqDotMinus_;
        
        //-
        volScalarField mLiqDotPlus_;


private:

    // Private Member Functions


        //- Disallow default bitwise copy construct
        schnerrSauerEvaporationSource
        (
            const schnerrSauerEvaporationSource&
        );

        //- Disallow default bitwise assignment
        void operator=(const schnerrSauerEvaporationSource&);
        
        //-
        dimensionedScalar alphaNuc(const dimensionedScalar& nNuc) const;
        
        //-
        tmp<volScalarField> rRb (const volScalarField& limAlphaLiq, const dimensionedScalar& nNuc);


public:

    //- Runtime type information
    TypeName("schnerrSauerEvaporationSource");


    // Constructors

        //- Construct from components
        schnerrSauerEvaporationSource
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~schnerrSauerEvaporationSource()
    {}


    // Member Functions

        virtual bool alwaysApply() const;


        // Public Functions

            //- Source term to fvMatrix<scalar>
            virtual void addSup(fvMatrix<scalar>& eqn, const label fieldI);
            
            //- Scalar
            virtual void addSup
            (
                const volScalarField& rho,
                fvMatrix<scalar>& eqn,
                const label fieldI
            );



        // I-O

            //- Write data
            virtual void writeData(Ostream&) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
