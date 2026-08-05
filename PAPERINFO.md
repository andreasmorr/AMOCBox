# PAPERINFO: AMOCBox

## Purpose

This file is the paper-facing summary for the AMOCBox submodule. Use it as the
local source of truth when updating `../main.tex`, especially the Methods,
Appendix, synthesis table, and figure captions. It intentionally omits general
usage instructions and keeps only information needed to describe the scientific
experiment.

If this file conflicts with `../main.tex`, prefer this file for the AMOCBox
experiment setup.

## Module Role In The Paper

AMOCBox provides the conceptual-model member of the AMOC resilience hierarchy. It
uses the reduced three-box AMOC model from Wood/Alkhayuon-style salinity-box
dynamics, represented here as a two-dimensional ODE in North Atlantic and
Tropical Atlantic salinity anomalies.

The submodule contributes:

- a CO2-parameter continuation of the AMOC on-state and other attractors,
- non-local resilience metrics around the AMOC-on attractor,
- the AMOCBox rows of the cross-model resilience synthesis,
- the AMOCBox paper figure showing representative trajectories and basin geometry.

## Source Files

The implemented experiment is defined by:

- `scripts/amoc3box_co2_continuation.jl`
- `scripts/amoc3box_export_paper_data.jl`
- `scripts/plotting_paper.py`
- `src/used_dynamical_systems.jl`
- `src/attractor_centered_stability.jl`

Primary paper-facing outputs are in:

- `data/paper/resilience_vs_co2_boxmodel.csv`
- `data/paper/basin_{co2}ppm.csv`
- `data/paper/trajectories_{co2}ppm.csv`
- `data/paper/attractors_{co2}ppm.csv`
- `plots/amocbox_paper.png`

## Model State And Units

The dynamical state is two-dimensional:

- `S_N`: North Atlantic salinity anomaly.
- `S_T`: Tropical Atlantic salinity anomaly.

Internal state units are `100 * (S - S_0)`, with `S_0 = 0.035` in the parameter
file convention. For paper plots, internal salinity `x` is displayed as
`35 + 10*x` psu. Therefore:

- `1.0` internal salinity unit equals `10 psu`.
- `0.1` internal salinity units equals `1 psu`.
- `0.2` internal salinity units equals `2 psu`.
- `0.005` internal salinity units equals `0.05 psu`.

The model also contains prescribed or diagnostic salinities for other boxes:

- Southern-box salinity `S_S` is prescribed by the parameter vector.
- Deep-box salinity `S_B` is prescribed by the parameter vector.
- Indo-Pacific salinity `S_IP` is diagnostically determined by salt conservation.

Do not describe the implemented AMOCBox reduction as one where the Southern
Ocean salinity is the diagnostic salt-conservation variable. In this code, the
diagnostic salinity is `S_IP`.

The AMOC strength diagnostic `q` is returned in cubic meters per second and is
converted to Sverdrups by dividing by `1e6`.

## Model Equations As Implemented

The ODE state is `S = (S_N, S_T)`. The freshwater hosing parameter `H` is present
in the model code but is set to `0.0` in both the 1xCO2 and 2xCO2 parameter
vectors used for the CO2 continuation.

Before evaluating the right-hand side, the hosing-adjusted fluxes are

```text
F_N <- F_N + 0.1311e6 * H
F_T <- F_T + 0.6961e6 * H
```

With `H = 0`, these reduce to the prescribed freshwater fluxes.

The prescribed dimensional salinities `S_S` and `S_B` are converted into the
internal anomaly convention:

```text
S_S <- 100 * (S_S - S_0)
S_B <- 100 * (S_B - S_0)
```

Indo-Pacific salinity is diagnosed from total salt content:

```text
S_IP = 100 * (
    C
    - (V_N*S_N + V_T*S_T + V_S*S_S + V_B*S_B) / 100
    - S_0 * (V_N + V_T + V_S + V_IP + V_B)
) / V_IP
```

The overturning strength is

```text
q = lambda * (alpha*(T_S - T_0) + beta*(S_N - S_S)/100)
    / (1 + lambda*alpha*mu)
```

The flow direction switches at `q = 0`.

For `q >= 0`:

```text
dS_N/dt = Y/V_N * (
    q*(S_T - S_N)/100
    + K_N*(S_T - S_N)/100
    - F_N*S_0
)

dS_T/dt = Y/V_T * (
    q*(gamma*S_S + (1 - gamma)*S_IP - S_T)/100
    + K_S*(S_S - S_T)/100
    + K_N*(S_N - S_T)/100
    - F_T*S_0
)
```

For `q < 0`:

```text
dS_N/dt = Y/V_N * (
    -q*(S_B - S_N)/100
    + K_N*(S_T - S_N)/100
    - F_N*S_0
)

dS_T/dt = Y/V_T * (
    -q*(S_N - S_T)/100
    + K_S*(S_S - S_T)/100
    + K_N*(S_N - S_T)/100
    - F_T*S_0
)
```

## Numerical Solver And Attractor Search

The dynamical system is constructed as `CoupledODEs(amoc_rule, u0, params)` with:

- solver: `Rosenbrock23()`,
- absolute tolerance: `1e-7`,
- relative tolerance: `1e-7`,
- initial condition: `SVector(0.0, 0.4)`.

The fixed attractor-search grid used for global continuation is:

- `S_N` in `[-1.0, 1.0]`, 101 grid points,
- `S_T` in `[-0.1, 2.0]`, 106 grid points.

The recurrence mapper options in the continuation script are:

- transient time `Ttr = 0.0`,
- sampling step `Delta t = 1.0`,
- `stop_at_Delta t = true`,
- `horizon_limit = 1e2`,
- `consecutive_lost_steps = 10000`,
- `consecutive_recurrences = 10000`,
- `consecutive_attractor_steps = 1000`,
- `consecutive_basin_steps = 1000`.

The proximity threshold used by the main resilience calculation is
`epsilon = 0.005` internal salinity units, i.e. `0.05 psu`.

## CO2 Forcing Coordinate

The experiment does not prescribe independent model configurations for every
CO2 level. It constructs a linear parameter curve between the 1xCO2 and 2xCO2
parameter vectors:

```text
params(t) = params_1x + t * (params_2x - params_1x)
```

The coordinate mapping used for labels is:

```text
CO2(t) = 280 ppm + t * (560 ppm - 280 ppm)
       = 280 ppm + 280 ppm * t
```

Key values:

- `t = 0`: 1xCO2, 280 ppm.
- `t = 1`: 2xCO2, 560 ppm.
- `t > 1`: linear extrapolation beyond the 2xCO2 parameter vector.

Nominal continuation settings:

- `T_START = 0.0`,
- `T_END = 2.4`,
- `T_STEP = 0.05`,
- 49 nominal continuation points from 280 ppm to 952 ppm.

In the current cached paper outputs:

- complete paper-facing resilience metrics are available from 280 ppm through
  896 ppm, corresponding to `t = 0.0` through `t = 2.2`;
- the AMOC-strength diagnostic in `resilience_vs_co2_boxmodel.csv` extends one
  additional step, through 910 ppm (`t = 2.25`);
- basin, trajectory, and attractor CSVs are exported from 280 ppm through
  896 ppm.

When writing the paper, do not state that all resilience metrics are available
through 952 ppm.

## CO2-Varying Parameters

The continuation linearly interpolates/extrapolates all parameters whose 1xCO2
and 2xCO2 values differ.

Changing parameters include:

- box volumes `V_N`, `V_T`, `V_S`, `V_IP`, `V_B`,
- freshwater fluxes `F_N`, `F_T`,
- temperatures `T_S`, `T_0`,
- mixing coefficients `K_N`, `K_S`,
- hydraulic constant `lambda`,
- mixing fraction `gamma`,
- momentum-advection coefficient `mu`,
- total salinity content `C`.

Parameters that are unchanged between the 1xCO2 and 2xCO2 vectors include:

- hosing `H = 0.0`,
- prescribed salinities `S_S`, initial `S_IP` placeholder, and `S_B`,
- expansion/contraction coefficients `alpha`, `beta`,
- reference salinity `S_0`,
- time-scaling parameter `Y`.

## Resilience Sampling

The continuation first identifies attractors using `global_continuation` with
`AttractorsViaRecurrences` and `AttractorSeedContinueMatch`.

The AMOC-on attractor ID is chosen as the attractor with the highest AMOC
strength at the first nonempty continuation step. Subsequent resilience metrics
are evaluated for this on-attractor while it remains present.

For the resilience metrics, the code does not sample from the fixed global
attractor-search grid. Instead, it uses a custom attractor-centered sampler:

- at each CO2 step, compute the centroid of the AMOC-on attractor;
- build a local rectangular sampling box centered on that centroid;
- half-width: `grid_delta = 0.2` internal salinity units, i.e. `+/-2 psu`;
- grid spacing used to define the sampling region: `grid_step = 0.02` internal
  salinity units, i.e. `0.2 psu`;
- random samples per parameter step: `RESILIENCE_SAMPLES = 10000`;
- sampling distribution: `Attractors.EverywhereUniform()`;
- finite-time horizon for stability measures: `FINITE_TIME = 1000.0`;
- distance object: `StrictlyMinimumDistance()`;
- proximity threshold: `epsilon = 0.005` internal salinity units, i.e. `0.05 psu`.

This setup means that the reported non-local resilience metrics describe the
AMOC-on basin in a local `+/-2 psu` salinity-perturbation box around the current
on-state, not in a fixed global rectangle.

## Paper-Facing Metrics

The synthesis CSV exports the following measures for the AMOC-on attractor:

- `amoc_strength_sv`: AMOC strength in Sv.
- `characteristic_return_time`: Attractors.jl characteristic return time.
- `local_resilience`: computed as `1 / characteristic_return_time`.
- `mean_convergence_time`: mean convergence time for sampled initial conditions.
- `basin_stability`: fraction of sampled initial conditions assigned to the
  AMOC-on attractor.
- `minimal_critical_shock_magnitude`: smallest sampled state-space distance from
  the on-attractor to an initial condition assigned away from the on-attractor.

Units and conversions:

- `amoc_strength_sv` is in Sverdrups.
- time metrics are in the model-year units used by the scripts.
- salinity-distance metrics are in internal salinity units unless converted;
  multiply by `10` to express them in psu.

For paper terminology:

- `basin_stability` is the implemented basin-volume/basin-fraction metric within
  the sampled `+/-2 psu` box.
- `minimal_critical_shock_magnitude` is the implemented minimal-critical-shock
  proxy from the finite sample.
- `local_resilience` is the inverse characteristic return time exported from the
  Attractors.jl stability measures. This is the intended local-resilience value
  for the AMOCBox synthesis.

## Paper Figure Data

The AMOCBox paper figure is produced by `scripts/plotting_paper.py`.

Current plotted scenarios:

- 280 ppm (`t = 0.0`),
- 448 ppm (`t = 0.6`).

For each exported CO2 level, `scripts/amoc3box_export_paper_data.jl` writes:

- `attractors_{co2}ppm.csv`: on/off attractor centroids and AMOC strength,
- `basin_{co2}ppm.csv`: a `60 x 60` basin grid,
- `trajectories_{co2}ppm.csv`: four representative perturbed trajectories.

The figure basin grid is centered on the on-attractor with:

- grid half-width `BASIN_DELTA = 0.2` internal salinity units (`+/-2 psu`),
- grid resolution `GRID_RES = 60` in each direction.

The four representative trajectories start near the four corners of this local
box:

- offset `0.8 * BASIN_DELTA = 0.16` internal salinity units,
- equivalent to `+/-1.6 psu` from the on-attractor center in each plotted
  salinity coordinate.

Trajectory output settings:

- integration length `T_MAX = 3000.0` model years,
- output interval `DT_TRAJ = 1.0` model year.

The attractor-neighborhood circles drawn in the plot have radius `0.005`
internal salinity units, i.e. `0.05 psu`.

Implementation note: the export helper constructs systems through
`amoc3box_at_t(t)`, whose proximity threshold is `epsilon = 0.01` internal units
for the exported basin classification. The main resilience metrics use
`epsilon = 0.005`.

## Recommended Paper Wording

Short methods wording:

> We analyse the reduced three-box AMOC model as a two-dimensional ODE in North
> Atlantic and Tropical Atlantic salinity anomalies. The model is continued along
> a linear parameter path between published 1xCO2 and 2xCO2 parameter vectors,
> labelled by the coordinate CO2(t) = 280 ppm + 280 ppm t. At each continuation
> step, we identify the AMOC-on attractor as the attractor with the largest
> overturning strength and evaluate resilience metrics from initial conditions
> sampled uniformly in a +/-2 psu box centered on that attractor.

Metric wording:

> Basin volume is represented by the fraction of sampled initial conditions in
> this on-state-centered perturbation box that return to the AMOC-on attractor.
> Mean convergence time and minimal critical shock are computed from the same
> sampled ensemble using the Attractors.jl stability-measure accumulator. Local
> resilience is reported as the inverse characteristic return time.

Figure-caption wording:

> AMOCBox trajectories and basin geometry are shown at 280 ppm and 448 ppm along
> the continuation. Phase-space axes are the North Atlantic and Tropical Atlantic
> salinity coordinates, converted from internal anomaly units to psu by
> `S_plot = 35 + 10 S_internal`.

## Do Not Say

- Do not say the AMOCBox resilience metrics were computed every 30 ppm or only
  over 285-480 ppm; the implemented continuation uses 14 ppm nominal spacing
  from the `t` step of 0.05.
- Do not say the complete resilience curves extend to 952 ppm; the current
  complete paper-facing metric range ends at 896 ppm.
- Do not describe the perturbation sampling for resilience metrics as a fixed
  global grid. It is attractor-centered at every CO2 step.
- Do not describe the implemented salt-conservation diagnostic as Southern Ocean
  salinity. In the code, `S_IP` is diagnostic while `S_S` and `S_B` are
  prescribed.
- Do not describe `local_resilience` as a separately hand-coded Jacobian
  calculation in this repository. It is exported as `1 / characteristic_return_time`.

