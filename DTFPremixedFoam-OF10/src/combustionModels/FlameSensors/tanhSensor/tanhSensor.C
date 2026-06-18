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

\*---------------------------------------------------------------------------*/

#include "tanhSensor.H"
#include "addToRunTimeSelectionTable.H"
#include "combustionModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace flameSensorModels
{
    defineTypeNameAndDebug(tanhSensor, 0);

    addToRunTimeSelectionTable
    (
        flameSensor,
        tanhSensor,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameSensorModels::tanhSensor::tanhSensor
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    flameSensor(dict, mesh),
    coeffsDict_(dict.subDict("tanhCoeffs")),
    beta_(coeffsDict_.lookup<scalar>("beta")),
    n_filters_(coeffsDict_.lookup<scalar>("n_filters")),
    sFilter_(mesh_),
    // HLG 2026 - Added: read optional qdotFloor guard parameters with safe defaults
    // so existing cases without these entries compile and run unchanged. Fix 6.
    qdotFloor_(coeffsDict_.lookupOrDefault<scalar>("qdotFloor", 0.0)),
    useQdotFloor_(coeffsDict_.lookupOrDefault<bool>("useQdotFloor", false))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::flameSensorModels::tanhSensor::correct()
{
    // Look up the combustion model by name to retrieve Qdot
    const word modelName
    (
        IOobject::groupName
        (
            combustionModel::combustionPropertiesName,
            ""
        )
    );

    tmp<volScalarField> tqdot =
        mesh_.lookupObject<combustionModel>(modelName).Qdot();

    // Apply simple box filter n_filters_ times to smooth Qdot
    for (int i = 0; i < n_filters_; i++)
    {
        tqdot = sFilter_(tqdot);
    }

    volScalarField qdot = Foam::mag(tqdot.ref());

    /*
    forAll(qdot, cellI)
    {
        qdot[cellI] = max(qdot[cellI], 100);
    }
    */

    const scalar qdotMax = gMax(qdot);

    Info << "  [tanhSensor] qdot : min = " << gMin(qdot)
         <<                    "  max = " << qdotMax << nl;

    // Inactive sensor when flame is absent: zero S_ and return early.
    // Two guards, in order of precedence:
    //   1. qdotFloor (optional, user-controlled): deactivate when gMax < floor.
    //   2. SMALL fallback (always active): prevents division by zero in cold flow.
    if (useQdotFloor_ && qdotMax < qdotFloor_)
    {
        S_ = dimensionedScalar(dimless, 0.0);
        Info << "  [tanhSensor] gMax(qdot) below floor " << qdotFloor_
             << " — sensor inactive" << nl;
        return;
    }

    if (qdotMax < SMALL)
    {
        S_ = dimensionedScalar(dimless, 0.0);
        return;
    }

    S_ = tanh(beta_ * qdot / dimensionedScalar(qdot.dimensions(), qdotMax));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flameSensorModels::tanhSensor::~tanhSensor()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// ************************************************************************* //

