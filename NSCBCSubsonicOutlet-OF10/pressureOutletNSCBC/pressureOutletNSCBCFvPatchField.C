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

// NSCBC subsonic outlet — pressure field
//
// Based on Poinsot & Lele (1992), J. Comput. Phys. 101, 104-129.
//
// At a subsonic outlet one wave enters the domain: the upstream-running
// acoustic wave L1 (wave speed u-c, Eq. 19). Its amplitude cannot be read
// from the interior and is instead estimated by relaxing the boundary
// pressure toward a target far-field value pInf (Section 3.2):
//
//     L1 = K*(p - pInf),   K = etaAc*(1 - M^2)*c/lInf          (*)
//
// The outgoing acoustic wave L5 (wave speed u+c, Eq. 23) is computed
// directly from one-sided interior gradients.
//
// The boundary pressure is advanced using the LODI pressure relation
// (Eq. 25):
//
//     dp/dt + (L1 + L5)/2 = 0
//
// The relaxation coefficient K in (*) follows the recommendation of Rudy &
// Strikwerda (1980): K ~ (1-M^2)*c/L, with etaAc as the tuning prefactor.
// A value etaAc = 0.25 is a widely used starting point.

#include "pressureOutletNSCBCFvPatchField.H"
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
Foam::pressureOutletNSCBCFvPatchField<Type>::pressureOutletNSCBCFvPatchField
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
Foam::pressureOutletNSCBCFvPatchField<Type>::pressureOutletNSCBCFvPatchField
(
    const pressureOutletNSCBCFvPatchField& ptf,
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
    pInf_(ptf.pInf_),
    gamma_(ptf.gamma_),
    etaAc_(ptf.etaAc_),
    fieldInf_(ptf.fieldInf_),
    lInf_(ptf.lInf_)
{}


template<class Type>
Foam::pressureOutletNSCBCFvPatchField<Type>::pressureOutletNSCBCFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    psiName_(dict.lookupOrDefault<word>("psi", "thermo:psi")),
    pInf_(readScalar(dict.lookup("pInf"))),
    gamma_(dict.lookupOrDefault<scalar>("gamma", 0.0)),
    etaAc_(readScalar(dict.lookup("etaAc"))),
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
Foam::pressureOutletNSCBCFvPatchField<Type>::pressureOutletNSCBCFvPatchField
(
    const pressureOutletNSCBCFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptpsf, iF),
    UName_(ptpsf.UName_),
    phiName_(ptpsf.phiName_),
    rhoName_(ptpsf.rhoName_),
    psiName_(ptpsf.psiName_),
    pInf_(ptpsf.pInf_),
    gamma_(ptpsf.gamma_),
    etaAc_(ptpsf.etaAc_),
    fieldInf_(ptpsf.fieldInf_),
    lInf_(ptpsf.lInf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::pressureOutletNSCBCFvPatchField<Type>::advectionSpeed() const
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
Foam::pressureOutletNSCBCFvPatchField<Type>::soundSpeed() const
{
    const fvPatchField<scalar>& psip =
        this->patch().template lookupPatchField<volScalarField, scalar>(psiName_);

    if (gamma_ > 0)
    {
        return sqrt(gamma_/psip);
    }

    if (this->db().template foundObject<fluidThermo>(physicalProperties::typeName))
    {
        const fluidThermo& thermo =
            this->db().template lookupObject<fluidThermo>(physicalProperties::typeName);
        const fvPatchScalarField& Tp =
            this->patch().template lookupPatchField<volScalarField, scalar>("T");
        return sqrt(thermo.gamma(Tp, this->patch().index())/psip);
    }

    FatalErrorInFunction
        << "gamma not specified on patch " << this->patch().name()
        << " and no fluidThermo found in registry.\n"
        << "Add 'gamma' to the boundary condition dictionary."
        << exit(FatalError);
    return tmp<scalarField>(nullptr);
}



template<class Type>
void Foam::pressureOutletNSCBCFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvMesh& mesh = this->internalField().mesh();

    word ddtScheme
    (
        mesh.schemes().ddt(this->internalField().name())
    );
    scalar deltaT = this->db().time().deltaTValue();

    const GeometricField<Type, fvPatchField, volMesh>& field =
        this->db().objectRegistry::template
        lookupObject<GeometricField<Type, fvPatchField, volMesh>>
        (
            this->internalField().name()
        );

    // Calculate the advection speed of the field wave
    // If the wave is incoming set the speed to 0.
    // const scalarField w(Foam::max(advectionSpeed(), scalar(0)));

    const  scalarField cP(soundSpeed());
    const  scalarField aP(advectionSpeed());

    const fvPatchScalarField& rhop =
            this->patch().template lookupPatchField<volScalarField, scalar>
            (
                rhoName_
            );

    const fvPatchVectorField& Up =
            this->patch().template lookupPatchField<volVectorField, vector>(UName_);



    label patchi = this->patch().index();

    // Non-reflecting outflow boundary
    // If lInf_ defined setup relaxation to the value fieldInf_.
    if (lInf_ > 0)
    {

	    if
		    (
		     ddtScheme == fv::EulerDdtScheme<scalar>::typeName
		     || ddtScheme == fv::CrankNicolsonDdtScheme<scalar>::typeName
		    )
		    {
			    // K = etaAc*(1-M^2)*c/lInf — relaxation rate (P&L Sec. 3.2)
			    const scalarField K(etaAc_*(1.0-sqr(aP/cP))*cP/lInf_);

			    // valueFraction balances the implicit L1 relaxation against
			    // the outgoing acoustic wave speed (u+c)/2 (P&L Eq. 25)
			    this->valueFraction() = (1.0 + K*deltaT/2.0)/(1.0 + K*deltaT/2.0 + (aP+cP)/2.0*deltaT*this->patch().deltaCoeffs()) ;
			    // refValue: p_old advanced by L1 relaxation source (P&L Eq. 25)
			    this->refValue() =
				    (
				     field.oldTime().boundaryField()[patchi]
				     + K * deltaT/2.0 * pInf_ * fieldInf_
				    )/( 1.0 + K * deltaT/2.0);
			    // refGrad: acoustic term from L5 = (u+c)*(dp/dn + rho*c*du/dn),
			    // projected onto the normal: contributes -rho*c*(du/dn) (P&L Eq. 23)
			    this->refGrad() = - 1.0 * rhop * cP * (this-> patch().nf() & Up.snGrad()) * fieldInf_;
		    }
	    else if (ddtScheme == fv::backwardDdtScheme<scalar>::typeName)
	    {
		    // K = etaAc*(1-M^2)*c/lInf — relaxation rate (P&L Sec. 3.2)
		    const scalarField K(etaAc_*(1.0-sqr(aP/cP))*cP/lInf_);

		    this->valueFraction() = (1.5 + K*deltaT/2.0)/(1.5 + K*deltaT/2.0 + (aP+cP)/2.0*deltaT*this->patch().deltaCoeffs()) ;
		    // refValue: BDF2 time level p_old advanced by L1 relaxation source
		    this->refValue() =
			    (
			     2.0*field.oldTime().boundaryField()[patchi]
                           - 0.5*field.oldTime().oldTime().boundaryField()[patchi]
			     + K * deltaT/2.0 * pInf_ * fieldInf_
			    )/( 1.5 + K * deltaT/2.0);
		    // refGrad: acoustic term from L5 (P&L Eq. 23)
		    this->refGrad() = - 1.0 * rhop * cP * (this-> patch().nf() & Up.snGrad()) * fieldInf_;
	    }
        else
        {
            FatalErrorInFunction
                    << "lInf_ must above 0 "
                    << exit(FatalError);
        }
    }

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::pressureOutletNSCBCFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    if (phiName_ != "phi")
    {
        os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
    }
    if (rhoName_ != "rho")
    {
        os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
    }
    if (UName_ != "U")
    {
        os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("pInf") << pInf_ << token::END_STATEMENT << nl;
    if (gamma_ > 0)
    {
        os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
    }
    os.writeKeyword("etaAc") << etaAc_ << token::END_STATEMENT << nl;

    if (lInf_ > 0)
    {
        os.writeKeyword("fieldInf") << fieldInf_ << token::END_STATEMENT << nl;
        os.writeKeyword("lInf") << lInf_ << token::END_STATEMENT << nl;
    }

    writeEntry(os, "value", *this);
}


// ************************************************************************* //
