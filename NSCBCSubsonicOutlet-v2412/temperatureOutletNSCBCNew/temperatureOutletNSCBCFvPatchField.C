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

// NSCBC subsonic outlet — temperature field (v2412 port).
//
// Wave amplitudes (P&L notation):
//
//   L1 = K*(p - pInf),   K = etaAc*(1-M^2)*c/lInf   [P&L Eq. 40]
//   L2 = u*(c^2*drho/dn - dp/dn)                     [P&L Eq. 20]
//   L5 = (u+c)*(dp/dn + rho*c*du/dn)                 [P&L Eq. 23]
//
// Temperature LODI [P&L Eq. 29]:
//
//     dT/dt + T/(rho*c^2)*[-(L2) + (gamma-1)/2*(L5+L1)] = 0
//
// Implementation note: using the ideal-gas relation, the entropy term
//     T/(rho*c^2)*L2 = -u*dT/dn + u*(gamma-1)/(gamma*rho*R)*dp/dn,
// so L2 is carried by the implicit convective term (valueFraction, speed u)
// together with the pressure-gradient refGrad. It is therefore NOT added
// explicitly to refValue; doing so would double-count the entropy wave.

#include "temperatureOutletNSCBCFvPatchField.H"
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
Foam::temperatureOutletNSCBCFvPatchField<Type>::temperatureOutletNSCBCFvPatchField
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
Foam::temperatureOutletNSCBCFvPatchField<Type>::temperatureOutletNSCBCFvPatchField
(
    const temperatureOutletNSCBCFvPatchField& ptf,
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
Foam::temperatureOutletNSCBCFvPatchField<Type>::temperatureOutletNSCBCFvPatchField
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
Foam::temperatureOutletNSCBCFvPatchField<Type>::temperatureOutletNSCBCFvPatchField
(
    const temperatureOutletNSCBCFvPatchField& ptpsf,
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
Foam::temperatureOutletNSCBCFvPatchField<Type>::advectionSpeed() const
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
Foam::temperatureOutletNSCBCFvPatchField<Type>::soundSpeed() const
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
void Foam::temperatureOutletNSCBCFvPatchField<Type>::updateCoeffs()
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

    const fvPatchField<scalar>& psip =
        this->patch().template lookupPatchField<volScalarField, scalar>(psiName_);

    const fvPatchScalarField& Tp =
        this->patch().template lookupPatchField<volScalarField, scalar>("T");

    const scalarField cP(soundSpeed());
    const scalarField aP(advectionSpeed());

    const scalarField gammaP(sqr(cP)*psip);
    const scalarField Rp(1.0/(psip*Tp));

    const fvPatchScalarField& rhop =
        this->patch().template lookupPatchField<volScalarField, scalar>(rhoName_);

    const fvPatchVectorField& Up =
        this->patch().template lookupPatchField<volVectorField, vector>(UName_);

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
            const scalarField L1(K*pp - K*pInf_);
            const scalarField L5
            (
                (aP+cP)*(pp.snGrad() + rhop*cP*(this->patch().nf() & Up.snGrad()))
            );

            // The entropy wave L2 = u*(c^2*drho/dn - dp/dn) is NOT added
            // explicitly here. Via the ideal-gas relation it is identically
            //     T/(rho c^2)*L2 = -u*dT/dn + u*(gamma-1)/(gamma*rho*R)*dp/dn,
            // i.e. it is already represented by the implicit convective term
            // (valueFraction below, advection speed u) plus the pressure-
            // gradient refGrad. Adding L2 again would double-count it.
            this->valueFraction() = 1.0/(1.0 + aP*deltaT*this->patch().deltaCoeffs());

            this->refValue() =
                field.oldTime().boundaryField()[patchi]
              - deltaT*(gammaP - 1.0)/gammaP*0.5*(L5+L1)/rhop/Rp*fieldInf_;

            this->refGrad() = (gammaP-1.0)/gammaP*pp.snGrad()/rhop/Rp*fieldInf_;
        }
        else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
        {
            const scalarField K(etaAc_*(1.0-sqr(aP/cP))*cP/lInf_);
            const scalarField L1(K*pp - K*pInf_);
            const scalarField L5
            (
                (aP+cP)*(pp.snGrad() + rhop*cP*(this->patch().nf() & Up.snGrad()))
            );

            // See note in the Euler branch: the entropy wave L2 is represented
            // by the implicit convective term + pressure-gradient refGrad and
            // must not be added explicitly (it would be double-counted).
            this->valueFraction() = 1.5/(1.5 + aP*deltaT*this->patch().deltaCoeffs());

            this->refValue() =
                (
                    2.0*field.oldTime().boundaryField()[patchi]
                  - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
                  - deltaT*(gammaP - 1.0)/gammaP*0.5*(L5+L1)/rhop/Rp*fieldInf_
                )/1.5;

            this->refGrad() = (gammaP-1.0)/gammaP*pp.snGrad()/rhop/Rp*fieldInf_;
        }
        else
        {
            FatalErrorInFunction
                << "unsupported ddt scheme '" << ddtScheme
                << "' in temperatureOutletNSCBC on patch " << this->patch().name()
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
void Foam::temperatureOutletNSCBCFvPatchField<Type>::write(Ostream& os) const
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
