# NSCBC Subsonic Outlet — OpenFOAM 10

Non-Reflecting Subsonic Characteristic Boundary Conditions (NSCBC) for compressible reacting flow in OpenFOAM 10. Three outlet boundary conditions are provided, one per transported variable.

## Boundary conditions

| Type | Applied to | Purpose |
|---|---|---|
| `pressureOutletNSCBC` | `p` | Non-reflecting pressure outlet |
| `temperatureOutletNSCBC` | `T` | Non-reflecting temperature outlet |
| `velocityOutletNSCBC` | `U` | Non-reflecting velocity outlet |

All three use the same relaxation approach: incoming acoustic wave amplitude is relaxed toward a target far-field state at a rate controlled by `etaAc` and `lInf`. Outgoing waves leave without reflection.

## Build

Compile each BC independently from its directory:

```bash
cd pressureOutletNSCBC  && wmake
cd temperatureOutletNSCBC && wmake
cd velocityOutletNSCBC  && wmake
```

Libraries land in `$FOAM_USER_LIBBIN`.

## Usage

Add the three libraries to `system/controlDict`:

```cpp
libs
(
    "libpressureOutletNSCBC.so"
    "libtemperatureOutletNSCBC.so"
    "libvelocityOutletNSCBC.so"
);
```

Then apply to the outlet patch in each field file. `gamma` is optional — if omitted it is read from the `fluidThermo` object at runtime, which is correct for mixtures where the composition varies across the patch.

### `0/p`
```cpp
outlet
{
    type        pressureOutletNSCBC;
    pInf        101325;     // target far-field pressure [Pa]
    etaAc       0.25;       // acoustic relaxation coefficient
    lInf        0.05;       // relaxation length scale [m]
    fieldInf    1;          // scaling factor on far-field term
    value       uniform 101325;
}
```

### `0/T`
```cpp
outlet
{
    type        temperatureOutletNSCBC;
    pInf        101325;
    etaAc       0.25;
    lInf        0.05;
    fieldInf    1;
    value       uniform 300;
}
```

### `0/U`
```cpp
outlet
{
    type        velocityOutletNSCBC;
    pInf        101325;
    etaAc       0.25;
    lInf        0.05;
    fieldInf    (1 0 0);    // outward unit normal of the outlet patch
    value       uniform (0 0 0);
}
```

If for any reason `fluidThermo` is not available (e.g. a non-reacting solver that does not register a thermo object), specify `gamma` explicitly:

```cpp
gamma   1.4;
```

## Parameters

| Parameter | Required | Description |
|---|---|---|
| `pInf` | yes | Far-field reference pressure [Pa] |
| `etaAc` | yes | Acoustic relaxation coefficient, typically 0.25 |
| `lInf` | yes | Relaxation length scale [m] |
| `fieldInf` | yes (when `lInf` > 0) | Direction multiplier on the characteristic terms. Scalar `1` for `p`/`T`; outward unit normal (e.g. `(1 0 0)`) for `U` |
| `gamma` | no | Ratio of specific heats. Auto-read from `fluidThermo` if omitted |
| `phi` | no | Flux field name (default: `phi`) |
| `rho` | no | Density field name (default: `rho`) |
| `U` | no | Velocity field name (default: `U`) |
| `p` | no | Pressure field name (default: `p`) |
| `psi` | no | Compressibility field name (default: `thermo:psi`) |

## Time schemes

`Euler`, `CrankNicolson` (first-order LODI integration), and `backward` (BDF2
time levels with the same spatial operators) are supported. Other `ddt` schemes
are rejected at run time.

## Documentation

See [`NSCBC_manual.pdf`](NSCBC_manual.pdf) for the full theory and the
term-by-term mapping of the implementation onto Poinsot & Lele (1992) and Baum,
Poinsot & Thévenin (1994), including how the LODI relations are realised as
OpenFOAM `mixed` (Robin) boundary conditions. The LaTeX source is
[`NSCBC_manual.tex`](NSCBC_manual.tex).

Species mass fractions at the outlet should use OpenFOAM's `advective`
boundary condition.
