/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::combustionModels::DTF

Description
    Dynamic Thickened Flame (DTF) LES combustion model.

    Implements the DTF approach of Legier, Poinsot & Veynante (2000), in which
    a local thickening factor F is active only near the flame front (F > 1) and
    relaxes to F = 1 away from it, preserving correct mixing outside the flame.

    The laminar flame speed (SL) and flame thickness (deltaL) are supplied
    via DTFCoeffs in combustionProperties. An efficiency function E accounts
    for unresolved sub-grid flame wrinkling.

    The effective reaction rate and heat release are scaled by E/F.

    Reference:
        Legier, J.P., Poinsot, T., Veynante, D. (2000).
        Dynamically thickened flame LES model for premixed and non-premixed
        turbulent combustion. CTR Summer Program 2000, pp. 157-168.

SourceFiles
    DTF.C

\*---------------------------------------------------------------------------*/

#ifndef DTF_H
#define DTF_H

#include "laminar.H"
#include "flameSensor.H"
#include "efficiencyFunction.H"
#include "autoPtr.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                            Class DTF Declaration
\*---------------------------------------------------------------------------*/

class DTF
:
    public laminar
{
    // Private Data
        scalar Fmax_;
        //- Target number of grid cells across the flame front (F = N*delta/deltaL)
        scalar N_;
        //- Sensor-modulated thickening factor (F=1 outside flame, F=Fs inside)
        volScalarField F_;
        //- Combined field written to disk: EF = E * F, used by transport models
        volScalarField EF_;
        //- Raw cell-size-based thickening factor (pre-sensor, capped at Fmax)
        volScalarField Fs_;
        //- Flame sensor: detects reaction zone, returns S in [0,1]
        autoPtr<flameSensor> flameSensor_;
        //- Efficiency function E: corrects for loss of flame-turbulence interaction due to thickening
        autoPtr<efficiencyFunction> efficiencyFunction_;
        //- Laminar flame thickness [m], read from DTFCoeffs
        scalar deltaL_;
        //- Laminar flame speed [m/s], read from DTFCoeffs
        scalar SL_;
    public:

    //- Runtime type information
    TypeName("DTF");


    // Constructors

        //- Construct from components
        DTF
        (
            const word& modelType,
            const fluidReactionThermo& thermo,
            const compressibleMomentumTransportModel& turb,
            const word& combustionProperties
        );

        //- Disallow default bitwise copy construction
        DTF(const DTF&);


        //- Destructor
        virtual ~DTF();


    // Member Functions

        //- Correct combustion rate
        virtual void correct();

        //- Fuel consumption rate matrix.
        virtual tmp<fvScalarMatrix> R(volScalarField& Y) const;

        //- Heat release rate [kg/m/s^3]
        virtual tmp<volScalarField> Qdot() const;

        scalar deltaL() const { return deltaL_; }
        scalar SL() const { return SL_; }

        const volScalarField& F() const { return F_; }
        

        //- Update properties from given dictionary
        virtual bool read();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const DTF&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
