/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

// NSCBC subsonic outlet — velocity field (v2412 port).
//
// Incoming acoustic wave [P&L Eq. 40]:   L1 = K*(p - pInf)
// Outgoing acoustic wave [P&L Eq. 23]:   L5 = (u+c)*(dp/dn + rho*c*du/dn)
//
// Velocity LODI — normal component [P&L Eq. 26]:
//
//     du/dt + (L5 - L1)/(2*rho*c) = 0
//
// Acoustic gradient relation for the refGrad [P&L Eq. 35]:
//
//     du/dn = -(1/(rho*c)) * dp/dn

#include "velocityOutletNSCBCFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "fluidThermo.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::velocityOutletNSCBCFvPatchField<Type>::velocityOutletNSCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    psiName_("thermo:psi"),
    pName_("p"),
    pInf_(0.0),
    gamma_(0.0),
    etaAc_(0.25),
    fieldInf_(Zero),
    lInf_(-GREAT)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


template<class Type>
Foam::velocityOutletNSCBCFvPatchField<Type>::velocityOutletNSCBCFvPatchField
(
    const velocityOutletNSCBCFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    psiName_(ptf.psiName_),
    pName_(ptf.pName_),
    pInf_(ptf.pInf_),
    gamma_(ptf.gamma_),
    etaAc_(ptf.etaAc_),
    fieldInf_(ptf.fieldInf_),
    lInf_(ptf.lInf_)
{}


template<class Type>
Foam::velocityOutletNSCBCFvPatchField<Type>::velocityOutletNSCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    UName_(dict.getOrDefault<word>("U", "U")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    psiName_(dict.getOrDefault<word>("psi", "thermo:psi")),
    pName_(dict.getOrDefault<word>("p", "p")),
    pInf_(dict.get<scalar>("pInf")),
    gamma_(dict.getOrDefault<scalar>("gamma", 0.0)),
    etaAc_(dict.get<scalar>("etaAc")),
    fieldInf_(Zero),
    lInf_(-GREAT)
{
    if (dict.found("value"))
    {
        fvPatchField<Type>::operator=
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fvPatchField<Type>::operator=(this->patchInternalField());
    }

    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;

    if (dict.readIfPresent("lInf", lInf_))
    {
        dict.lookup("fieldInf") >> fieldInf_;

        if (lInf_ < 0.0)
        {
            FatalIOErrorInFunction(dict)
                << "unphysical lInf specified (lInf < 0)" << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << exit(FatalIOError);
        }
    }
}


template<class Type>
Foam::velocityOutletNSCBCFvPatchField<Type>::velocityOutletNSCBCFvPatchField
(
    const velocityOutletNSCBCFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptpsf, iF),
    UName_(ptpsf.UName_),
    phiName_(ptpsf.phiName_),
    rhoName_(ptpsf.rhoName_),
    psiName_(ptpsf.psiName_),
    pName_(ptpsf.pName_),
    pInf_(ptpsf.pInf_),
    gamma_(ptpsf.gamma_),
    etaAc_(ptpsf.etaAc_),
    fieldInf_(ptpsf.fieldInf_),
    lInf_(ptpsf.lInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::velocityOutletNSCBCFvPatchField<Type>::advectionSpeed() const
{
    const surfaceScalarField& phi =
        this->db().objectRegistry::template lookupObject<surfaceScalarField>
        (phiName_);

    const fvsPatchField<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField, scalar>
        (
            phiName_
        );

    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchScalarField& rhop =
            this->patch().template lookupPatchField<volScalarField, scalar>
            (
                rhoName_
            );

        return phip/(rhop*this->patch().magSf());
    }
    else
    {
        return phip/this->patch().magSf();
    }
}

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::velocityOutletNSCBCFvPatchField<Type>::soundSpeed() const
{
    const fvPatchField<scalar>& psip =
        this->patch().template lookupPatchField<volScalarField, scalar>(psiName_);

    if (gamma_ > 0)
    {
        return sqrt(gamma_/psip);
    }

    if (this->db().template foundObject<fluidThermo>(fluidThermo::typeName))
    {
        const fluidThermo& thermo =
            this->db().template lookupObject<fluidThermo>(fluidThermo::typeName);
        const scalarField gammaP
        (
            thermo.gamma()().boundaryField()[this->patch().index()]
        );
        return sqrt(gammaP/psip);
    }

    FatalErrorInFunction
        << "gamma not specified on patch " << this->patch().name()
        << " and no fluidThermo found in registry.\n"
        << "Add 'gamma' to the boundary condition dictionary."
        << exit(FatalError);
    return tmp<scalarField>(nullptr);
}


template<class Type>
void Foam::velocityOutletNSCBCFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();

    word ddtScheme
    (
        mesh.ddtScheme(this->internalField().name())
    );
    scalar deltaT = this->db().time().deltaTValue();

    const GeometricField<Type, fvPatchField, volMesh>& field =
        this->db().objectRegistry::template
        lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            this->internalField().name()
        );

    const  scalarField cP(soundSpeed());
    const  scalarField aP(advectionSpeed());

    const fvPatchScalarField& rhop =
        this->patch().template lookupPatchField<volScalarField, scalar>
        (
            rhoName_
        );

    const fvPatchScalarField& pp =
        this->patch().template lookupPatchField<volScalarField, scalar>(pName_);

    label patchi = this->patch().index();

    if (lInf_ > 0)
    {
        if
        (
            ddtScheme == fv::EulerDdtScheme<scalar>::typeName
         || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
        )
        {
            const scalarField K(etaAc_*(1.0-sqr(aP/cP))*cP/lInf_);

            this->valueFraction() = 1.0/(1.0 + (aP+cP)/2.0*deltaT*this->patch().deltaCoeffs());
            this->refValue() =
                (
                    field.oldTime().boundaryField()[patchi]
                  + K * deltaT/2.0 /rhop / cP  * (pp - pInf_) * fieldInf_
                );
            this->refGrad() = - 1.0 / rhop / cP * pp.snGrad() * fieldInf_;
        }
        else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
        {
            const scalarField K(etaAc_*(1.0-sqr(aP/cP))*cP/lInf_);

            this->valueFraction() = 1.5/(1.5 + (aP+cP)/2.0*deltaT*this->patch().deltaCoeffs());
            this->refValue() =
                (
                    2.0*field.oldTime().boundaryField()[patchi]
                  - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
                  + K * deltaT/2.0 / rhop / cP * (pp - pInf_) * fieldInf_
                )/1.5;
            this->refGrad() = - 1.0 / rhop / cP * pp.snGrad() * fieldInf_;
        }
        else
        {
            FatalErrorInFunction
                << "unsupported ddt scheme '" << ddtScheme
                << "' in velocityOutletNSCBC on patch " << this->patch().name()
                << exit(FatalError);
        }
    }
    else
    {
        FatalErrorInFunction
            << "lInf_ must be > 0 on patch " << this->patch().name()
            << exit(FatalError);
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::velocityOutletNSCBCFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    if (phiName_ != "phi")
    {
        os.writeEntry("phi", phiName_);
    }
    if (rhoName_ != "rho")
    {
        os.writeEntry("rho", rhoName_);
    }
    if (UName_ != "U")
    {
        os.writeEntry("U", UName_);
    }
    if (pName_ != "p")
    {
        os.writeEntry("p", pName_);
    }
    os.writeEntry("pInf", pInf_);
    if (gamma_ > 0)
    {
        os.writeEntry("gamma", gamma_);
    }
    os.writeEntry("etaAc", etaAc_);

    if (lInf_ > 0)
    {
        os.writeEntry("fieldInf", fieldInf_);
        os.writeEntry("lInf", lInf_);
    }

    this->writeEntry("value", os);
}


// ************************************************************************* //
