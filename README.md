# AMOCBox

Resilience analysis of the Alkhayuon et al. 3-box AMOC model under increasing CO₂, using [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/).

The model is parametrised for 1×CO₂ and 2×CO₂ settings. Interpolating between these two parameter vectors gives a continuum of AMOC states under increasing CO₂. Stability and resilience measures are computed along this continuum.

---

## File structure

```
AMOCBox/
├── scripts/
│   ├── amoc3box_co2_continuation.jl   # Global continuation along CO₂ parameter curve
│   ├── amoc3box_export_paper_data.jl  # Export CSVs for paper figures
│   └── plotting_paper.py              # Paper figure (reads exported CSVs)
├── src/
│   └── used_dynamical_systems.jl      # Box model definitions (1×CO₂, 2×CO₂ variants)
├── data/
│   ├── co2_continuation/              # Cached results from CO₂ continuation runs (.jld2)
│   └── paper/                         # CSV exports for plotting
│       ├── basin_{1xco2,2xco2}.csv        # 60×60 basin label grid
│       ├── trajectories_{1xco2,2xco2}.csv # Perturbed IC trajectories
│       ├── attractors_{1xco2,2xco2}.csv   # Attractor positions in salinity space
│       └── resilience_vs_co2_boxmodel.csv # Resilience measures vs CO₂
├── plots/
│   └── amocbox_paper.pdf              # Output paper figure
├── Project.toml
└── Manifest.toml
```

---

## Scripts

### `amoc3box_co2_continuation.jl`

Runs global attractor continuation along `params(t) = params_1x + t * (params_2x - params_1x)`, where `t=0` is pre-industrial and `t=1` is 2×CO₂. Stability and resilience measures are computed at each step and saved via DrWatson's `produce_or_load`. Also exports `data/paper/resilience_vs_co2_boxmodel.csv` with AMOC strength (Sv) and resilience metrics along the CO₂ parameter curve.

### `amoc3box_export_paper_data.jl`

Exports all data required by `plotting_paper.py` to `data/paper/`. Computes a 60×60 basin grid via `AttractorsViaRecurrences`, selects three representative initial conditions (on-basin, boundary-crossing, off-basin) per CO₂ level, integrates trajectories for 3 000 time units, and writes basin grids, trajectories, and attractor positions as CSVs.

### `src/used_dynamical_systems.jl`

Box model definitions for the 1×CO₂ and 2×CO₂ parameter sets, including the `amoc_strength` diagnostic (returns overturning in m³/s; divide by 1×10⁶ for Sv).

---

## Paper figure

`plotting_paper.py` produces a publication-quality 4-panel figure using the shared design language from `../amoc_plot_style.py`.

**Figure layout:**
- **Top row** (shorter): AMOC strength vs time for 3 selected trajectories at 1×CO₂ (left) and 2×CO₂ (right). ON/OFF equilibria shown as dashed horizontal lines.
- **Bottom row** (square): 2D salinity phase portrait (N. Atlantic vs Tropical Atlantic box salinity). Basin of attraction shown as pale background colour. Trajectories are time-shaded (alpha increases with time). Attractor positions marked with stars.

Output: `plots/amocbox_paper.pdf`

---

## Usage

From the project root:

```bash
# 1. Run continuation analysis
julia --project scripts/amoc3box_co2_continuation.jl

# 2. Export paper data (basin grids, trajectories, attractors)
julia --project scripts/amoc3box_export_paper_data.jl

# 3. Generate paper figure
python scripts/plotting_paper.py
```

Results are cached in `data/` and figures are written to `plots/`.

---

## Dependencies

Julia with [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/), [DynamicalSystems.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/dynamicalsystems/stable/), [Attractors.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/attractors/stable/), and [CairoMakie](https://makie.org/). Python requires `numpy` and `matplotlib`. Install Julia packages with:

```julia
using Pkg; Pkg.instantiate()
```
