/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
 
#include "DTF.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(DTF, 0);
    addToRunTimeSelectionTable(combustionModel, DTF, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::DTF::DTF
(
    const word& modelType,
    const fluidReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    laminar(modelType, thermo, turb, combustionProperties),
    Fmax_(this->coeffs().template lookup<scalar>("Fmax")),
    N_(this->coeffs().template lookup<scalar>("N")),
    F_
     (
         IOobject
         (
             "F",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
         ),
         this->mesh(),
         dimensionedScalar(dimless, 1)
     ),
    EF_
     (
         IOobject
         (
             "EF",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
         ),
         this->mesh(),
         dimensionedScalar(dimless, 1)
     ),
    Fs_
     (
         IOobject
         (
             "Fs",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::AUTO_WRITE
         ),
         this->mesh(),
         dimensionedScalar(dimless, 1)
     ),
    flameSensor_(flameSensor::New(coeffs(),turb.mesh())),
    efficiencyFunction_(efficiencyFunction::New(turb.mesh(),coeffs(),turb,Fmax_, *this)),
    deltaL_(this->coeffs().template lookup<scalar>("deltaL")),
    SL_(this->coeffs().template lookup<scalar>("SL"))
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::DTF::~DTF()
{}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::DTF::correct()
{
    Info << nl << "// * * * * * * * * DTF::correct() * * * * * * * * //" << endl;

    laminar::correct();
    flameSensor_->correct();

    // Compute cell-size-based thickening factor Fs and sensor-blended F BEFORE
    // efficiencyFunction_->correct(), so the efficiency function (e.g.
    // durandPolifke) sees the current step's F via dtfModel_.F() rather than
    // the previous step's value.  Without this ordering, EF_ would mix a
    // stale E with a fresh F.
    // delta_g = V^(1/3): cube root cell size, consistent with cubeRootVolDelta
    // and correct for 1D/2D/3D unstructured meshes (nGeometricD() is wrong for 1D/2D).
    forAll(F_, celli)
    {
        const scalar delta_g = pow(mesh_.V()[celli], 1.0/3.0);
        Fs_[celli] = min(N_ * delta_g / deltaL_, Fmax_);
        F_[celli]  = max(1.0 + (Fs_[celli] - 1.0)*flameSensor_->S()[celli], 1.0);
    }
    Info << "Thickening factor Fs (raw, pre-sensor) : min = " << Foam::gMin(Fs_) << "  max = " << Foam::gMax(Fs_) << endl;
    Info << "Thickening factor F  (sensor-modulated) : min = " << Foam::gMin(F_) << "  max = " << Foam::gMax(F_) << endl;

    efficiencyFunction_->correct();

    EF_ = efficiencyFunction_->E() * F_;
    Info << "Efficiency function E : min = " << Foam::gMin(efficiencyFunction_->E()) << "  max = " << Foam::gMax(efficiencyFunction_->E()) << endl;

    Info << "// * * * * * * * * * * * * * * * * * * * * * * * * //" << nl << endl;
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::DTF::R(volScalarField& Y) const
{
    return efficiencyFunction_->E()/F_*laminar::R(Y);
}

Foam::tmp<Foam::volScalarField>
Foam::combustionModels::DTF::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + ":Qdot"),
        efficiencyFunction_->E()/F_*laminar::Qdot()
    );
}


bool Foam::combustionModels::DTF::read()
{
    if (laminar::read())
    {
        this->coeffs().lookup("Fmax") >> Fmax_;
        this->coeffs().lookup("N")    >> N_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
