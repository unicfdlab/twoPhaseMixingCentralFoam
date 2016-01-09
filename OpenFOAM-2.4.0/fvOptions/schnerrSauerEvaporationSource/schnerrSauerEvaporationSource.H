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
    Heat exchanger source model, in which the heat exchanger is defined as a
    selection of cells.

    The total heat exchange source is given by:
    \f[
        Q_t = e(\phi, \dot{m}_2) (T_2 - T_1) \phi c_p
    \f]

    where:
    \vartable
        Q_t       | total heat source
        e(\phi,\dot{m}_2) | effectivenes table
        \phi      | net mass flux entering heat exchanger [kg/s]
        \dot{m}_2 | secondary mass flow rate [kg/s]
        T_1       | primary inlet temperature [K]
        T_2       | secondary inlet temperature [K]
        c_p       | specific heat capacity [J/kg/K]
    \endvartable


    The distribution inside the hear exchanger is given by:
    \f[
        Q_c = \frac{V_c |U_c| (T_c - T_{ref})}{\sum(V_c |U_c| (T_c - T_{ref}))}
    \f]

    where:
    \vartable
        Q_c     | source for cell
        V_c     | volume of the cell [m3]
        U_c     | local cell velocity [m/s]
        T_c     | local call temperature [K]
        T_{ref} | min or max(T) in cell zone depending on the sign of Q_t [K]
    \endvartable

    \heading Source usage

    Example usage:
    \verbatim
    schnerrSauerEvaporationSource1
    {
        type            schnerrSauerEvaporationSource;
        active          true;
        selectionMode   cellZone;
        cellZone        porosity;

        schnerrSauerEvaporationSourceCoeffs
        {
            secondaryMassFlowRate   1.0;
            secondaryInletT         336;
            primaryInletT           293;
            faceZone                facesZoneInletOriented;
            outOfBounds             clamp;
            fileName                "effTable";
        }
    }
    \endverbatim

    The effectiveness table is described in terms of the primary and secondary
    mass flow rates.  For example, the table:

                           secondary MFR
                       |  0.1   0.2   0.3
                  -----+-----------------
                  0.02 |   A     B     C
     primary MFR  0.04 |   D     E     F
                  0.06 |   G     H     I


    Is specified by the following:

        (
            0.02
            (
                (0.1    A)
                (0.2    B)
                (0.3    C)
            ),
            0.04
            (
                (0.1    D)
                (0.2    E)
                (0.3    F)
            ),
            0.06
            (
                (0.1    G)
                (0.2    H)
                (0.3    I)
            )
        );


Note
- the table with name "fileName" should have the same units as the
  secondary mass flow rate and kg/s for phi
- faceZone is the faces at the inlet of the cellzone, it needs to be
  created with flip map flags. It is used to integrate the net mass flow
  rate into the heat exchanger


SourceFiles
    schnerrSauerEvaporationSource.C

\*---------------------------------------------------------------------------*/

#ifndef schnerrSauerEvaporationSource_H
#define schnerrSauerEvaporationSource_H

#include "fvOption.H"
#include "autoPtr.H"
#include "interpolation2DTable.H"

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
        dimensionedScalar n_;
        
        //- Nucleation site diameter
        dimensionedScalar dNuc_;
        
        //- Condensation constant
        dimensionedScalar Cc_;
        
        //- Vaporization constant
        dimensionedScalar Cv_;
	
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
        dimensionedScalar alphaNuc() const;
        
        //-
        tmp<volScalarField> rRb (const volScalarField& limAlphaLiq);


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
