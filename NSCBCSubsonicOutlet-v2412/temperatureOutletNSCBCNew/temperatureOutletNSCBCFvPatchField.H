/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::temperatureOutletNSCBCFvPatchField

Group
    grpOutletBoundaryConditions

Description
    Non-reflecting subsonic outlet for temperature, based on the Navier-Stokes
    Characteristic Boundary Conditions (NSCBC) of Poinsot & Lele (1992),
    J. Comput. Phys. 101, 104-129.

    Temperature is advanced using the LODI relation (Eq. 29):

        dT/dt + (T/rho/c^2) * [-(L2) + (gamma-1)/2*(L5 + L1)] = 0

    where:
      L1 = (u-c)*(dp/dn - rho*c*du/dn)  incoming acoustic  [relaxed to pInf]
      L2 = u*(c^2*drho/dn - dp/dn)       entropy wave       [from interior]
      L5 = (u+c)*(dp/dn + rho*c*du/dn)  outgoing acoustic  [from interior]

    The entropy wave L2 allows density/temperature variations to pass through
    the boundary without spurious reflection. The factor T/(rho*c^2) is
    evaluated as 1/(gamma*rho*R) using local thermo properties, so spatially
    varying heat capacities in reacting or multicomponent flows are handled
    correctly.

    Supported time schemes: Euler, CrankNicolson, backward.

Usage
    \table
        Property  | Description                              | Required | Default
        pInf      | far-field (target) pressure [Pa]         | yes      |
        etaAc     | acoustic relaxation coefficient (sigma)  | yes      |
        lInf      | relaxation length scale [m]              | yes      |
        fieldInf  | multiplier for BC direction (1 or -1)   | yes      |
        gamma     | heat capacity ratio (overrides thermo)   | no       |
        phi       | flux field name                          | no       | phi
        rho       | density field name                       | no       | rho
        p         | pressure field name                      | no       | p
        U         | velocity field name                      | no       | U
    \endtable

    Example:
    \verbatim
    outlet
    {
        type        temperatureOutletNSCBC;
        pInf        101325;
        etaAc       0.25;
        lInf        0.05;
        fieldInf    1;
        value       uniform 300;
    }
    \endverbatim

SourceFiles
    temperatureOutletNSCBCFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef temperatureOutletNSCBCFvPatchField_H
#define temperatureOutletNSCBCFvPatchField_H

#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class temperatureOutletNSCBCFvPatchField Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class temperatureOutletNSCBCFvPatchField
:
    public mixedFvPatchField<Type>
{
protected:

    // Private data
        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;


        //- Name of the density field used to normalise the mass flux
        //- if necessary
        word rhoName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Name of the pressure field
        word pName_;

        //- Far-field (target) pressure [Pa]
        scalar pInf_;

        //- Heat capacity ratio gamma = Cp/Cv (overrides thermo if > 0)
        scalar gamma_;

        //- Acoustic relaxation coefficient sigma (P&L Sec. 3.2)
        scalar etaAc_;

        //- Field value of the far-field
        Type fieldInf_;

        //- Relaxation length-scale
        scalar lInf_;

public:

    //- Runtime type information
    TypeName("temperatureOutletNSCBC");


    // Constructors

        //- Construct from patch and internal field
        temperatureOutletNSCBCFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        temperatureOutletNSCBCFvPatchField
        (
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given temperatureOutletNSCBCFvPatchField
        //  onto a new patch
        temperatureOutletNSCBCFvPatchField
        (
            const temperatureOutletNSCBCFvPatchField<Type>&,
            const fvPatch&,
            const DimensionedField<Type, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchField<Type>> clone
        (
            const DimensionedField<Type, volMesh>& iF
        ) const
        {
            return tmp<fvPatchField<Type>>
            (
                new temperatureOutletNSCBCFvPatchField<Type>(*this, iF)
            );
        }

        //- Construct as copy setting internal field reference
        temperatureOutletNSCBCFvPatchField
        (
            const temperatureOutletNSCBCFvPatchField&,
            const DimensionedField<Type, volMesh>&
        );


    // Member functions

        // Access

            //- Return the velocity at infinity
            const scalar& pInf() const
            {
                return pInf_;
            }

            //- Return reference to the velocity at infinity to allow adjustment
            scalar& pInf()
            {
                return pInf_;
            }

            //- Return the field at infinity
            const Type& fieldInf() const
            {
                return fieldInf_;
            }

            //- Return reference to the field at infinity to allow adjustment
            Type& fieldInf()
            {
                return fieldInf_;
            }

            //- Return the relaxation length-scale
            scalar lInf() const
            {
                return lInf_;
            }

            //- Return reference to the relaxation length-scale
            //  to allow adjustment
            scalar& lInf()
            {
                return lInf_;
            }

            //- Return the heat capacity ratio
            scalar gamma() const
            {
                return gamma_;
            }

            //- Return reference to the heat capacity ratio to allow adjustment
            scalar& gamma()
            {
                return gamma_;
            }

            //- Return the etaAc NSCBC ratio
            scalar etaAc() const
            {
                return etaAc_;
            }

            //- Return reference to the etaAc NSCBC ratio
            scalar& etaAc()
            {
                return etaAc_;
            }



        // Evaluation functions

            //- Calculate and return the advection speed at the boundary
            virtual tmp<scalarField> advectionSpeed() const;
            
            //- Calculate soundSpeed
            virtual tmp<scalarField> soundSpeed() const;

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "temperatureOutletNSCBCFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
