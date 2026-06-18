# NSCBC Subsonic Outlet — OpenFOAM v2412 (port)

Partially reflecting (weakly non-reflecting) subsonic outlet for compressible
reacting flow, ported from the OpenFOAM 10 version
([`../NSCBCSubsonicOutlet-OF10`](../NSCBCSubsonicOutlet-OF10)) to OpenFOAM
v2412. The physics and the boundary-condition formulation are identical; only
the library interfaces differ. Three outlet conditions are provided, one per
transported variable.

## Boundary conditions

| Type | Applied to | Purpose |
|---|---|---|
| `pressureOutletNSCBC` | `p` | Non-reflecting pressure outlet (relaxation to `pInf`) |
| `temperatureOutletNSCBC` | `T` | Non-reflecting temperature outlet |
| `velocityOutletNSCBC` | `U` | Non-reflecting velocity outlet |

The incoming acoustic wave amplitude is relaxed toward a target far-field
pressure at a rate set by `etaAc` and `lInf`; outgoing waves leave without
reflection. `gamma`, `c` and `R` are evaluated face-by-face from the
`fluidThermo` object, so the condition is correct across the strong property
gradients of a flame.

## Build

```bash
./Allwmake     # builds libpressureOutletNSCBCNew.so, etc. into $FOAM_USER_LIBBIN
```

## Usage

Register the libraries in `system/controlDict`:

```cpp
libs
(
    "libpressureOutletNSCBCNew.so"
    "libtemperatureOutletNSCBCNew.so"
    "libvelocityOutletNSCBCNew.so"
);
```

Apply all three on the outlet patch (with `advective` for each species `Yi`):

```cpp
// 0/p
outlet { type pressureOutletNSCBC;    pInf 101325; etaAc 0.25; lInf 0.05; fieldInf 1;       value uniform 101325; }
// 0/T
outlet { type temperatureOutletNSCBC; pInf 101325; etaAc 0.25; lInf 0.05; fieldInf 1;       value uniform 300; }
// 0/U
outlet { type velocityOutletNSCBC;    pInf 101325; etaAc 0.25; lInf 0.05; fieldInf (1 0 0); value uniform (0 0 0); }
```

`fieldInf` is the scalar `1` for `p`/`T` and the outward unit normal (e.g.
`(1 0 0)`) for `U`. Specify `gamma` explicitly only if no `fluidThermo` object
is registered.

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `pInf` | yes | Far-field reference pressure [Pa] |
| `etaAc` | yes | Acoustic relaxation coefficient, typically 0.25 |
| `lInf` | yes | Relaxation length scale [m] |
| `fieldInf` | yes (when `lInf` > 0) | Direction multiplier (scalar `1` for `p`/`T`; outward unit normal for `U`) |
| `gamma` | no | Ratio of specific heats. Auto-read from `fluidThermo` if omitted |
| `phi`, `rho`, `U`, `p`, `psi` | no | Field-name overrides |

## Time schemes

`Euler`, `CrankNicolson`, and `backward` (BDF2) are supported.

## Documentation

See [`NSCBC_manual.pdf`](NSCBC_manual.pdf) (source
[`NSCBC_manual.tex`](NSCBC_manual.tex)) for the full theory and the term-by-term
mapping onto Poinsot & Lele (1992) and Baum, Poinsot & Thévenin (1994).
