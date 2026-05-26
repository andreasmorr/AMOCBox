# AMOCBox

Resilience analysis of the Alkhayuon et al. 3-box AMOC model under increasing CO₂, using [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/).

The model is parametrized for 1×CO₂ and 2×CO₂ settings. Interpolating between (and extrapolating beyond) these two parameter vectors gives a continuum of AMOC states under increasing CO₂. Multiple stability and resilience measures are computed along this continuum. In a parallel experiment, freshwater hosing (H) is used as the continuation parameter at fixed 1×CO₂ and 2×CO₂ to map out the bifurcation structure.

## File structure

```
AMOCBox/
├── scripts/
│   ├── amoc3box_co2_continuation.jl   # Global continuation along CO2 parameter curve
│   └── amoc3box_hosing_scan.jl        # Hosing continuation at 1xCO2 and 2xCO2
├── src/
│   └── used_dynamical_systems.jl      # Box model definitions (1xCO2, 2xCO2 variants)
├── data/
│   ├── co2_continuation/              # Cached results from CO2 continuation runs (.jld2)
│   └── hosing_continuation/           # Cached results from hosing scan runs (.jld2)
├── Project.toml
└── Manifest.toml
```

## Scripts

### `amoc3box_co2_continuation.jl`
Runs global attractor continuation along the parameter curve
`params(t) = params_1x + t * (params_2x - params_1x)`, where `t=0` is pre-industrial and `t=1` is 2×CO₂. Stability and resilience measures are computed at each step and saved via DrWatson's `produce_or_load`.

### `amoc3box_hosing_scan.jl`
Runs global attractor continuation over freshwater hosing H = 0 → 0.55 Sv separately for the 1×CO₂ and 2×CO₂ parameter sets, producing a comparison of resilience measures across the two CO₂ scenarios.

## Usage

From the project root:

```bash
julia --project scripts/amoc3box_co2_continuation.jl
julia --project scripts/amoc3box_hosing_scan.jl
```

Results are cached in `data/` and figures are written to `plots/`.

## Dependencies

Julia with [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/), [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/), [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/), and [CairoMakie](https://makie.org/). Install with:

```julia
using Pkg; Pkg.instantiate()
```
