#ifndef durandPolifke_H
#define durandPolifke_H

#include "efficiencyFunction.H"
#include "simpleFilter.H"
#include "momentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace efficiencyFunctionModels
{

/*---------------------------------------------------------------------------*\
                    Class durandPolifke Declaration

    Colin (2000) efficiency function for the Dynamic Thickened Flame (DTF)
    model, using the TF' finite-volume SGS velocity operator of Durand &
    Polifke (2007) for compatibility with unstructured tetrahedral meshes.

    References:
      [1] Colin et al., Phys. Fluids 12(7):1843-1863, 2000.
      [2] Durand & Polifke, ASME GT2007-28188, 2007.

    ---- Efficiency function ------------------------------------------------

      E     = Xi(r0) / Xi(r1)                                    [2] Eq. 11
      Xi    = 1 + alpha * Gamma(r, s) * s                        [2] Eq.  7
      Gamma = 0.75 * exp(-1.2 / s^0.3) * r^(2/3)                [2] Eq. 10
      s     = u' / SL
      r0    = Delta_e / deltaL          (unthickened flame ratio)
      r1    = Delta_e / (F * deltaL)    (thickened flame ratio)

    ---- Test filter scale Delta_e ------------------------------------------

      [2] Eq. 9 recommends:
        Delta_e = delta_l^1 = F * deltaL
      so that the test filter sits exactly at the thickened flame scale.
      This is what the code implements.

      This differs from the common shortcut Delta_e = delta (LES grid width),
      which is simpler but does not match the paper's stated choice.

      With Delta_e = F*deltaL: r0 = F, r1 = 1 in every cell.

    ---- SGS velocity: TF' finite-volume operator ---------------------------

      u' = c * dx * |curl(U) - curl(U_filtered)|                 [2] Eq. 20

      dx = delta = V^(1/3)  (cubeRootVol local grid scale)

      U_filtered is one pass of the volume-weighted box filter:
        U_filtered_i = sum(U_k * V_k) / sum(V_k)  (face neighbours k)
                                                                  [2] Eq. 21

      OpenFOAM's simpleFilter performs this average.  On non-uniform meshes
      the volume-gradient cross-terms of [2] Eq. 23 are neglected — Durand
      notes this is exact on uniform meshes and acceptable elsewhere.

      This replaces the structured finite-difference operator ([1] Eq. 52 /
      [2] Eq. 16) which requires indirect Cartesian neighbours and cannot be
      applied to unstructured or tetrahedral meshes.

      Similarity constant: c ~ 2.0 (Cook 1997, cited in [2]).

    ---- Constant alpha -----------------------------------------------------

      [2] Eq. 8 gives: alpha = 2*ln2 / (3 * cms * (Re_t^0.5 - 1)), cms=0.28
      In practice alpha is read from the dictionary as a calibration constant
      because Re_t is not well-defined locally in a general LES.
      Typical range: alpha ~ 0.3-0.6 for Re_t = 100-10000.

    ---- Numerical guards (not in the paper) --------------------------------

      s < SMALL  ->  E = 1     (no turbulence, no correction needed)
      F < 1      ->  warning, F clamped to 1
      E clipped to [1, Emax],  Emax = Fmax^(2/3): the high-Re asymptote
        as s->inf, Xi(r0)/Xi(r1) -> r0^(2/3) = F^(2/3), so Fmax^(2/3) is
        the tightest physically consistent upper bound.

    ---- Dictionary --------------------------------------------------------

      durandPolifkeCoeffs
      {
          alpha   0.5;    // wrinkling constant; ~ 0.3-0.6, see [2] Eq. 8
          // c    2.0;    // TF' similarity constant (default 2.0, Cook 1997)
      }

\*---------------------------------------------------------------------------*/

class durandPolifke
:
    public efficiencyFunction
{

    // Private data

        dictionary coeffsDict_;

        //- Wrinkling model constant ([2] Eq. 8)
        scalar alpha_;

        //- TF' scale similarity constant ([2] Eq. 20, Cook 1997)
        scalar c_;

        //- Upper clip on E = Fmax^(2/3), the high-Re asymptote
        scalar Emax_;

        //- Counter for periodic diagnostic output
        label outputCounter_;

        const Foam::combustionModels::DTF& dtfModel_;

        //- SGS velocity magnitude [m/s] (AUTO_WRITE for diagnostics)
        volScalarField uPrime_;

        //- Finite-volume box filter (one pass = TF' test filter, [2] Eq. 21)
        simpleFilter sFilter_;


public:

    TypeName("durandPolifke");

    durandPolifke
    (
        const dictionary&,
        const fvMesh&,
        const compressibleMomentumTransportModel& turb,
        scalar F,
        const Foam::combustionModels::DTF& dtfModel
    );

    durandPolifke(const durandPolifke&) = delete;

    virtual ~durandPolifke();

    virtual void correct();

    void operator=(const durandPolifke&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace efficiencyFunctionModels
} // End namespace Foam

#endif

// ************************************************************************* //
