/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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
    Foam::turbulenceThermophysicalTransportModels::DTFunityLewisEddyDiffusivity

Description
    Eddy-diffusivity based energy gradient heat flux model for RAS or LES
    of turbulent flow.  Specie fluxes are computed assuming a unity turbulent
    Lewis number.

Usage
    \verbatim
    LES
    {
        model           DTFunityLewisEddyDiffusivity;
        Prt             0.85;
    }
    \endverbatim

SourceFiles
    DTFunityLewisEddyDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef DTFunityLewisEddyDiffusivity_H
#define DTFunityLewisEddyDiffusivity_H
#include "unityLewisEddyDiffusivity.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                       Class DTFunityLewisEddyDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceThermophysicalTransportModel>
class DTFunityLewisEddyDiffusivity
:
    public unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
{

const fvMesh& mesh_;

public:

    typedef typename TurbulenceThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        TurbulenceThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename TurbulenceThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("DTFunityLewisEddyDiffusivity");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        DTFunityLewisEddyDiffusivity
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        //- Construct from a type name, a momentum transport model and a thermo
        //  model, and whether to default the turbulent Prandtl number to one
        //  if it is not specified
        DTFunityLewisEddyDiffusivity
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo,
            const bool allowDefaultPrt = false
        );


    //- Destructor
    virtual ~DTFunityLewisEddyDiffusivity()
    {}


    // Member Functions

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
            const volScalarField& S  = mesh_.lookupObject<volScalarField>("S");

            // Molecular kappa is enhanced by EF (thickening * efficiency).
            // Turbulent alphat is suppressed inside the flame by (1-S): EF already
            // accounts for sub-grid flame-turbulence interaction via E, so retaining
            // full alphat in the flame zone would double-count sub-grid diffusion.
            return this->thermo().kappa()*EF
                 + (1-S)*this->alphat()*this->thermo().Cp();
        }


        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
            const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
            const volScalarField& S  = mesh_.lookupObject<volScalarField>("S");

            return this->thermo().kappa()/this->thermo().Cp()*EF
                 + (1-S)*this->alphat();
        }


        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const
        {
            const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
            const volScalarField& S  = mesh_.lookupObject<volScalarField>("S");

            return volScalarField::New
            (
                "DEff",
                this->thermo().kappa()/this->thermo().Cp()*EF + (1-S)*this->alphat()
            );
        }

        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
            const volScalarField& S  = mesh_.lookupObject<volScalarField>("S");
            return
                this->thermo().kappa().boundaryField()[patchi]
               *EF.boundaryField()[patchi]
              + (1-S.boundaryField()[patchi])
               *this->alphat(patchi)*this->thermo().Cp().boundaryField()[patchi];
        }

        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
            const volScalarField& S  = mesh_.lookupObject<volScalarField>("S");
            return
                this->thermo().kappa().boundaryField()[patchi]
               /this->thermo().Cp().boundaryField()[patchi]
               *EF.boundaryField()[patchi]
              + (1-S.boundaryField()[patchi])*this->alphat(patchi);
        }

        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const
        {
            const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
            const volScalarField& S  = mesh_.lookupObject<volScalarField>("S");
            return
                this->thermo().kappa().boundaryField()[patchi]
               /this->thermo().Cp().boundaryField()[patchi]
               *EF.boundaryField()[patchi]
              + (1-S.boundaryField()[patchi])*this->alphat(patchi);
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DTFunityLewisEddyDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
