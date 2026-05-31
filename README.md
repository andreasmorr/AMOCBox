# AMOCBox

Resilience analysis of the Alkhayuon et al. 3-box AMOC model under increasing CO₂, using [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/).

The model is parametrised for 1×CO₂ (280 ppm, pre-industrial) and 2×CO₂ (560 ppm) settings. The parameter curve is extrapolated beyond 2×CO₂ up to t=2.4 (≈952 ppm). Stability and resilience measures are computed along this continuum. The export script produces one set of CSVs per t step (49 steps total), each tagged by the rounded CO₂ value (e.g. `280ppm`, `308ppm`, …, `952ppm`). The paper figure compares two scenarios: 280 ppm (pre-industrial, t=0) and 896 ppm (t=2.2).

---

## File structure

```
AMOCBox/
├── scripts/
│   ├── amoc3box_co2_continuation.jl   # Global continuation along CO₂ parameter curve
│   ├── amoc3box_export_paper_data.jl  # Export CSVs for paper figures
│   └── plotting_paper.py              # Paper figure (reads exported CSVs)
├── src/
│   ├── used_dynamical_systems.jl      # Box model definitions (1×CO₂, 2×CO₂, 896 ppm variants)
│   └── attractor_centered_stability.jl # Custom stability measures with per-step centred sampling
├── data/
│   ├── co2_continuation/              # Cached results from CO₂ continuation runs (.jld2)
│   └── paper/                         # CSV exports for plotting
│       ├── basin_{co2}ppm.csv              # 60×60 basin label grid (one per t step)
│       ├── trajectories_{co2}ppm.csv       # Perturbed IC trajectories (one per t step)
│       ├── attractors_{co2}ppm.csv         # Attractor positions in salinity space (one per t step)
│       └── resilience_vs_co2_boxmodel.csv  # Resilience measures vs CO₂
├── plots/
│   └── amocbox_paper.png              # Output paper figure (300 dpi PNG)
├── Project.toml
└── Manifest.toml
```

---

## Scripts

### `amoc3box_co2_continuation.jl`

Runs global attractor continuation along `params(t) = params_1x + t * (params_2x - params_1x)`, where `t=0` is pre-industrial (280 ppm) and `t=1` is 2×CO₂ (560 ppm). The curve is extrapolated to `t=2.4` (≈952 ppm). Stability and resilience measures are computed at each step and saved via DrWatson's `produce_or_load`. Also directly exports `data/paper/resilience_vs_co2_boxmodel.csv` with AMOC strength (Sv) and the four resilience metrics used in the synthesis figure.

**Attractor-centred sampling**: the resilience measures use a custom function (`src/attractor_centered_stability.jl`) rather than `Attractors.stability_measures_along_continuation`. At each continuation step the centroid of the AMOC-on attractor is computed, and initial conditions for the stability measures are drawn uniformly from a ±2 psu box around that centroid (±0.2 in model units, where 1 model unit = 10 psu). This ensures the perturbation experiment is always physically anchored to the current on-state position as it drifts with CO₂, rather than sampling from a fixed global grid.

**Resilience measures exported to CSV**: `characteristic_return_time`, `mean_convergence_time`, `basin_stability`, `minimal_critical_shock_magnitude`.

### `amoc3box_export_paper_data.jl`

Loads the JLD2 continuation cache and exports one set of CSVs per t step to `data/paper/`, skipping steps where the on-attractor is absent. For each exported step (tagged `{co2}ppm`) it:
- uses `AttractorsViaProximity` with the pre-found continuation attractors (fast, no fresh attractor search),
- computes a 60×60 basin grid centred on the on-attractor (±0.2 model units = ±2 psu),
- integrates four trajectories from the grid corners (±0.8×0.2 model units),
- writes `basin_{co2}ppm.csv`, `trajectories_{co2}ppm.csv`, and `attractors_{co2}ppm.csv`.

Must be run after `amoc3box_co2_continuation.jl` has produced the JLD2 cache.

### `src/used_dynamical_systems.jl`

Box model definitions for the 1×CO₂, 2×CO₂, and 896 ppm parameter sets, including the `amoc_strength` diagnostic (returns overturning in m³/s; divide by 1×10⁶ for Sv). The 896 ppm parameter set uses temperatures extrapolated along the same 1x→2x direction at t=2.2.

---

## Paper figure

`plotting_paper.py` produces a publication-quality 4-panel figure using the shared design language from `../amoc_plot_style.py`.

**Figure layout:**
- **Top row** (shorter): AMOC strength vs time for 4 selected trajectories at 280 ppm (pre-industrial, left) and 896 ppm (right). ON/OFF equilibria shown as dashed horizontal lines.
- **Bottom row** (square): 2D salinity phase portrait (N. Atlantic vs Tropical Atlantic box salinity). Basin of attraction shown as pale background colour. Trajectories are time-shaded (alpha increases with time). Attractor positions marked with stars. Axes are labelled in psu: the model's internal salinity variable has offset 35 psu (i.e. a model value of 0 corresponds to 35 psu) and a scale of 10 psu per unit (i.e. a step of 0.1 in model units equals 1 psu).

Output: `plots/amocbox_paper.png` (300 dpi PNG)

---

## Usage

From the project root:

```bash
# 1. Run continuation analysis (also exports resilience_vs_co2_boxmodel.csv)
julia --project scripts/amoc3box_co2_continuation.jl

# 2. Export paper data (basin grids, trajectories, attractors for all t steps)
julia --project scripts/amoc3box_export_paper_data.jl

# 3. Generate paper figure
python scripts/plotting_paper.py
```

Results are cached in `data/` and figures are written to `plots/`.

---

## Dependencies

Julia with [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/), [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/), [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/), and [CairoMakie](https://makie.org/). Python requires `numpy`, `matplotlib`, and `pandas`. Install Julia packages with:

```julia
using Pkg; Pkg.instantiate()
```
