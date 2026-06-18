/*---------------------------------------------------------------------------*\
    Colin (2000) efficiency function with the TF' finite-volume SGS velocity
    operator of Durand & Polifke (2007) for unstructured tetrahedral meshes.

    References:
      [1] Colin et al., Phys. Fluids 12(7):1843-1863, 2000.
      [2] Durand & Polifke, ASME GT2007-28188, 2007.

    Key implementation choices and their paper basis:

      Delta_e = F * deltaL                    [2] Eq. 9  (paper's own choice)
      u' = c * dx * |curl(U) - curl(Uf)|     [2] Eq. 20 (TF' operator)
      Uf = volume-weighted box filter         [2] Eq. 21 (simpleFilter)
      E  = Xi(r0) / Xi(r1)                   [2] Eq. 11
      Xi = 1 + alpha * Gamma * s             [2] Eq.  7
      Gamma = 0.75*exp(-1.2/s^0.3)*r^(2/3)  [2] Eq. 10

    What is not exact vs the paper:

      - The volume-gradient cross-terms of [2] Eq. 23 are dropped (simpleFilter
        approximation). Durand notes this is exact on uniform meshes.
      - alpha is read from the dictionary rather than computed from Re_t via
        [2] Eq. 8. This is standard engineering practice.
      - F is a local sensor-modulated field, not a global constant. This is
        more correct for DTF but differs from the paper's simple presentation.

\*---------------------------------------------------------------------------*/

#include "durandPolifke.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcCurl.H"
#include "DTF.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace efficiencyFunctionModels
{
    defineTypeNameAndDebug(durandPolifke, 0);

    addToRunTimeSelectionTable
    (
        efficiencyFunction,
        durandPolifke,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::efficiencyFunctionModels::durandPolifke::durandPolifke
(
    const dictionary& dict,
    const fvMesh& mesh,
    const compressibleMomentumTransportModel& turb,
    scalar Fmax,
    const Foam::combustionModels::DTF& dtfModel
)
:
    efficiencyFunction(dict, mesh, turb, Fmax, dtfModel),
    coeffsDict_(dict.subDict("durandPolifkeCoeffs")),
    alpha_(coeffsDict_.lookup<scalar>("alpha")),
    c_(coeffsDict_.lookupOrDefault<scalar>("c", 2.0)),
    Emax_(Foam::pow(Fmax, 2.0/3.0)),
    outputCounter_(0),
    dtfModel_(dtfModel),
    uPrime_
    (
        IOobject
        (
            "durandPolifke_uPrime",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("durandPolifke_uPrime", dimLength/dimTime, 0.0)
    ),
    sFilter_(mesh_)
{
    Info<< "Durand-Polifke (TF') efficiency function initialised:" << nl
        << "  alpha    = " << alpha_ << nl
        << "  c        = " << c_     << nl
        << "  Emax     = " << Emax_  << " (= Fmax^(2/3), [2] p.3)" << nl
        << "  Delta_e  = F * deltaL             ([2] Eq. 9)" << nl
        << "  u' op    = TF' finite-volume       ([2] Eq. 20)" << nl
        << "  Suitable for unstructured / tetrahedral meshes." << nl
        << endl;

    if (alpha_ < 0.1)
    {
        WarningInFunction
            << "durandPolifke: alpha = " << alpha_ << " is unusually small." << nl
            << "  [2] Eq. 8 gives alpha ~ 0.3-0.6 for Re_t = 100-10000." << nl
            << "  This will produce a very weak efficiency correction." << nl
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::efficiencyFunctionModels::durandPolifke::~durandPolifke()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::efficiencyFunctionModels::durandPolifke::correct()
{
    const scalar deltaL = dtfModel_.deltaL();
    const scalar SL     = dtfModel_.SL();

    if (deltaL <= 0)
        FatalErrorInFunction
            << "durandPolifke: deltaL = " << deltaL << " <= 0." << exit(FatalError);

    if (SL <= 0)
        FatalErrorInFunction
            << "durandPolifke: SL = " << SL << " <= 0." << exit(FatalError);

    // Local grid scale dx = V^(1/3).  In OpenFOAM LES this is the cubeRootVol
    // delta field.  This is also the Δx in [2] Eq. 20.
    const volScalarField& dx =
        mesh_.lookupObject<volScalarField>("delta");

    const volScalarField& Fc = dtfModel_.F();

    // -----------------------------------------------------------------------
    // SGS velocity — TF' finite-volume operator
    // [2] Eq. 20:  u' = c * Δx * |curl(U) - curl(U_filtered)|
    //
    // U_filtered is the volume-weighted box filter over face-sharing
    // neighbours ([2] Eq. 21):
    //   U_filtered_i = sum(U_k * V_k) / sum(V_k)
    //
    // OpenFOAM's simpleFilter implements exactly this average.  On non-uniform
    // meshes the volume-gradient cross-terms of [2] Eq. 23 are neglected;
    // this is Durand's own stated simplification for practical use.
    //
    // Dimensional check:
    //   curl(U) - curl(U_filtered)  -> [1/s]
    //   dx                          -> [m]
    //   c * dx * |...|              -> [m/s]   correct
    // -----------------------------------------------------------------------

    volVectorField U_filtered
    (
        IOobject
        (
            "durandPolifke_U_filtered",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        sFilter_(turb_.U())   // one pass = TF' box filter ([2] Eq. 21)
    );

    uPrime_ = c_ * dx * mag(fvc::curl(turb_.U()) - fvc::curl(U_filtered));

    // -----------------------------------------------------------------------
    // Efficiency function
    // [2] Eq. 7, 10, 11
    //
    //   s     = u' / SL
    //   Gamma = 0.75 * exp(-1.2 / s^0.3) * r^(2/3)
    //   Xi    = 1 + alpha * Gamma * s
    //   E     = Xi(r0) / Xi(r1),   clipped to [1, Emax]
    //
    // Test filter scale: Delta_e = F * deltaL  ([2] Eq. 9)
    //   This is the paper's stated choice — it places the test filter at the
    //   thickened flame scale so the wrinkling ratio is physically meaningful.
    //
    //   r0 = Delta_e / deltaL        = F         (real flame ratio)
    //   r1 = Delta_e / (F * deltaL) = 1          (thickened flame ratio)
    //
    // Note: with this Delta_e choice, r1 = 1 always, so Xi(r1) is the same
    // in every cell (it only varies through s, which is cell-local).
    // Xi(r0) varies because F is local (sensor-modulated).
    //
    // E >= 1 analytically: r0 = F >= 1 = r1, so Xi(r0) >= Xi(r1).
    // The max(E, 1) guard is kept for floating-point safety only.
    // -----------------------------------------------------------------------

    bool fBelowOne = false;

    forAll(E_, celli)
    {
        const scalar Floc = Fc[celli];

        if (Floc < 1.0 - SMALL)
        {
            fBelowOne = true;
        }

        const scalar FlocSafe = max(Floc, 1.0);
        const scalar s        = uPrime_[celli] / SL;

        // No turbulence: no wrinkling correction needed
        if (s < SMALL)
        {
            E_[celli] = 1.0;
            continue;
        }

        // Delta_e = F * deltaL  ([2] Eq. 9)
        // r0 = Delta_e / deltaL      = F     (unthickened flame)
        // r1 = Delta_e / (F*deltaL)  = 1     (thickened flame, used implicitly
        //                                     in G1 below since r1^(2/3) = 1)
        const scalar r0 = FlocSafe;   // = (F*deltaL) / deltaL

        // Gamma ([2] Eq. 10) — exponential part shared between r0 and r1
        const scalar expFac = 0.75 * Foam::exp(-1.2 / Foam::pow(s, 0.3));

        const scalar G0  = expFac * Foam::pow(r0, 2.0/3.0);
        const scalar G1  = expFac; // r1^(2/3) = 1^(2/3) = 1

        const scalar Xi0 = 1.0 + alpha_ * G0 * s;
        const scalar Xi1 = 1.0 + alpha_ * G1 * s;

        // [2] Eq. 11
        E_[celli] = min(max(Xi0 / Xi1, 1.0), Emax_);
    }

    if (fBelowOne)
    {
        WarningInFunction
            << "durandPolifke: F < 1 detected in some cells." << nl
            << "  F clamped to 1 in efficiency calculation." << nl
            << "  Check your thickening factor field." << nl
            << endl;
    }

    if (outputCounter_++ % 10 == 0)
    {
        Info<< "Durand-Polifke (TF') efficiency function diagnostics:" << nl
            << "  uPrime : min = " << gMin(uPrime_)
            << "  max = "          << gMax(uPrime_) << " [m/s]" << nl
            << "  E      : min = " << gMin(E_)
            << "  max = "          << gMax(E_)      << " [-]"   << nl
            << endl;
    }
}


// ************************************************************************* //
