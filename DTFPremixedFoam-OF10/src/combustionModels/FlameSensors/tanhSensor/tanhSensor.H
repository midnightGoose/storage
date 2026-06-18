/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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
    Foam::flameSensorModels::tanhSensor

Description
    Hyperbolic tangent flame sensor.

    Detects the reaction zone by computing S = tanh(beta * |Qdot| / max|Qdot|),
    where Qdot is the heat release rate field filtered n_filters times with a
    simple box filter. The resulting sensor field S lies in [0, 1].

    Coefficients are read from the tanhCoeffs sub-dictionary:
    - beta      : sharpness parameter (higher = sharper transition)
    - n_filters : number of box filter passes applied to Qdot before sensing

SourceFiles
    tanhSensor.C

\*---------------------------------------------------------------------------*/

#ifndef tanhSensor_H
#define tanhSensor_H

#include "flameSensor.H"
#include "simpleFilter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flameSensorModels
{

/*---------------------------------------------------------------------------*\
                        Class tanhSensor Declaration
\*---------------------------------------------------------------------------*/

class tanhSensor
:
    public flameSensor
{
    // Private Data

        //- Sub-dictionary of tanhSensor coefficients
        dictionary coeffsDict_;

        //- Sharpness parameter for the tanh function
        scalar beta_;

        //- Number of box filter passes applied to Qdot
        scalar n_filters_;

        //- Simple box filter instance
        simpleFilter sFilter_;

        // HLG 2026 - Added: numerical stability guard for cold-flow / pre-ignition.
        // When gMax(qdot) is near zero the tanh normalisation is ill-conditioned.
        // qdotFloor_ is an absolute lower bound (e.g. 1e3 W/m3); if useQdotFloor_
        // is true and gMax(qdot) falls below it the sensor is zeroed for that step.
        // This is intentionally NOT a physics-based laminar reference — a fixed
        // Cantera-derived value would suppress the sensor in turbulent LES where
        // filtered qdot on a coarse grid is far below the laminar peak.
        //.
        scalar qdotFloor_;

        //- Switch to enable/disable the qdotFloor guard
        bool useQdotFloor_;


public:

    //- Runtime type information
    TypeName("tanh");


    // Constructors

        //- Construct from dictionary and mesh
        tanhSensor
        (
            const dictionary&,
            const fvMesh&
        );

        //- Disallow default bitwise copy construction
        tanhSensor(const tanhSensor&) = delete;


    //- Destructor
    virtual ~tanhSensor();

    //- Update the sensor field S_
    void correct() override;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const tanhSensor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace flameSensorModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

