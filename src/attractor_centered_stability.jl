"""
    attractor_centered_stability.jl

Custom replacement for `Attractors.stability_measures_along_continuation` that
samples initial conditions from a grid **centred on the on-state attractor** at
each parameter continuation step, rather than from a fixed global grid.

Motivation: the global grid used by `AttractorsViaRecurrences` to find attractors
covers a large region of state space. For resilience measures we want to evaluate
how robustly the system stays in the on-state when perturbed from near that state.
Using a locally centred grid of ±`grid_delta` (model units) around the on-attractor
makes the sampling physically meaningful and consistent across CO₂ levels, even as
the attractor drifts in state space.

Unit convention (3-box AMOC model):
  model unit = 100 × (S − S₀),  S₀ = 0.035 psu (35 psu)
  ⟹  1 model unit = 10 psu,  0.1 model unit = 1 psu
  ⟹  ±2 psu  ≡  ±0.2 model units  (default `grid_delta`)
"""

using Statistics: mean
using Attractors: ProgressMeter

"""
    stability_measures_attractor_centered(
        ds, attractors_cont, pcurve, on_id;
        grid_delta, grid_step,
        ε, weighting_distribution, finite_time, distance,
        samples_per_parameter, proximity_mapper_options,
        show_progress
    ) → Dict{String, Vector{Dict{Int,Float64}}}

Compute stability / resilience measures along a parameter continuation curve,
sampling initial conditions from a box centred on the on-state attractor at each
step.

# Arguments
- `ds`: the `DynamicalSystem` (parameters are mutated in-place per step).
- `attractors_cont`: vector of attractor dicts from `global_continuation`.
- `pcurve`: vector of parameter-update lists, one per continuation step.
- `on_id`: attractor ID that corresponds to the AMOC-on state (from `find_on_id`).

# Keyword arguments
- `grid_delta = 0.2`: half-width of the sampling box in model units (0.2 ≡ ±2 psu).
- `grid_step = 0.02`: grid spacing in model units (same resolution as the global grid).
- `ε`, `weighting_distribution`, `finite_time`, `distance`: identical semantics to
  `stability_measures_along_continuation` — each may be a scalar, an `AbstractVector`
  (indexed by step), or a `Function(p, attractors)`.
- `samples_per_parameter`: number of ICs sampled per continuation step.
- `proximity_mapper_options`: keyword args forwarded to `AttractorsViaProximity`.
- `show_progress`: show a ProgressMeter bar.

# Returns
Same format as `stability_measures_along_continuation`:
`Dict` mapping measure name → `Vector{Dict{attractor_id → value}}`.
Steps where the on-attractor is absent return empty inner dicts.
"""
function stability_measures_attractor_centered(
    ds,
    attractors_cont,
    pcurve,
    on_id::Int;
    grid_delta::Float64 = 0.2,
    grid_step::Float64  = 0.02,
    ε = nothing,
    weighting_distribution = Attractors.EverywhereUniform(),
    finite_time = 1.0,
    distance    = Attractors.StrictlyMinimumDistance(),
    samples_per_parameter::Int  = 1000,
    proximity_mapper_options    = NamedTuple(),
    show_progress::Bool = true,
)
    n        = length(pcurve)
    progress = ProgressMeter.Progress(
        n; desc = "Attractor-centred stability measures:", enabled = show_progress,
    )

    measures_cont = []
    measure_names = nothing

    for (i, p) in enumerate(pcurve)
        set_parameters!(ds, p)
        attractors = attractors_cont[i]

        # Skip steps where the on-attractor has disappeared
        if isempty(attractors) || !haskey(attractors, on_id) || isempty(attractors[on_id])
            push!(measures_cont, Dict{String, Dict{Int64, Float64}}())
            ProgressMeter.next!(progress)
            continue
        end

        # ── Centroid of on-state attractor ────────────────────────────────────
        A_on    = attractors[on_id]
        centroid = vec(mean(Matrix(A_on); dims = 1))
        sn_c, st_c = centroid[1], centroid[2]

        # ── Build ±grid_delta sampling grid centred on the on-attractor ──────
        xg = range(sn_c - grid_delta, sn_c + grid_delta; step = grid_step)
        yg = range(st_c - grid_delta, st_c + grid_delta; step = grid_step)
        sampler, = statespace_sampler((xg, yg))

        # ── Per-step parameter dispatch (scalar / vector / function) ─────────
        ε_ = _resolve(ε,                    i, p, attractors)
        wd = _resolve(weighting_distribution, i, p, attractors)
        ft = _resolve(finite_time,            i, p, attractors)
        d  = _resolve(distance,               i, p, attractors)

        # ── Accumulate stability measures ─────────────────────────────────────
        accumulator = StabilityMeasuresAccumulator(
            AttractorsViaProximity(ds, attractors; ε = ε_, proximity_mapper_options...);
            weighting_distribution = wd, finite_time = ft, distance = d,
        )

        basins_fractions(accumulator, sampler; N = samples_per_parameter, show_progress = false)
        measures, _ = finalize_accumulator(accumulator)

        if measure_names === nothing
            measure_names = collect(keys(measures))
        end
        push!(measures_cont, measures)
        ProgressMeter.next!(progress)
    end

    # ── Transpose to expected output format ───────────────────────────────────
    transposed = Dict{String, Vector{Dict{Int64, Float64}}}()
    isnothing(measure_names) && return transposed
    for name in measure_names
        transposed[name] = Vector{Dict{Int64, Float64}}()
    end
    for measures in measures_cont
        for name in measure_names
            push!(transposed[name], get(measures, name, Dict{Int64, Float64}()))
        end
    end
    return transposed
end

# ─── Internal helper ─────────────────────────────────────────────────────────

"""Resolve a per-step parameter that may be a scalar, vector, or function."""
function _resolve(x, i::Int, p, attractors)
    x isa AbstractVector && return x[i]
    x isa Function       && return x(p, attractors)
    return x
end
