"""
    amoc3box_co2_continuation.jl

Global continuation of the 3-box AMOC model along a linear parameter curve
that interpolates from the 1xCO2 parameter vector to the 2xCO2 vector and
then extrapolates further in the same direction.

Parameter curve:
  params(t) = params_1x + t * (params_2x - params_1x)
  t = 0   → 1xCO2 (pre-industrial baseline)
  t = 1   → 2xCO2 (~2°C warmer)
  t > 1   → extrapolation beyond 2xCO2

Attractors are tracked by global continuation; stability / resilience measures
are computed at each step and plotted in a single multi-panel figure.

Run from the project root:
    julia --project scripts/amoc3box_co2_continuation.jl
"""

using DrWatson
@quickactivate "AMOCResilience"

include(srcdir("used_dynamical_systems.jl"))
include(srcdir("attractor_centered_stability.jl"))

using DynamicalSystems
using Attractors
using CairoMakie
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# Tunable parameters
# ─────────────────────────────────────────────────────────────────────────────

const T_START = 0.0    # 1xCO2
const T_END   = 2.4    # beyond 2xCO2
const T_STEP  = 0.05

const CONTINUATION_SAMPLES = 100
const RESILIENCE_SAMPLES   = 1_000
const FINITE_TIME          = 1000.0

const COMPUTE_STABILITY = true   # set false to skip stability measures (faster)
const PLOT_BASINS       = false  # set true to compute and save per-step basin plots

const MEASURE_YLABELS = Dict(
    "basin_fraction"                        => "Basin fraction",
    "basin_stability"                       => "Basin stability",
    "finite_time_basin_stability"           => "Finite-time basin stab.",
    "mean_convergence_time"                 => "Mean conv. time (yr)",
    "maximal_convergence_time"              => "Max conv. time (yr)",
    "median_convergence_time"               => "Median conv. time (yr)",
    "mean_convergence_pace"                 => "Mean conv. pace",
    "maximal_convergence_pace"              => "Max conv. pace",
    "median_convergence_pace"               => "Median conv. pace",
    "minimal_critical_shock_magnitude"      => "Min. critical shock",
    "mean_noncritical_shock_magnitude"      => "Mean noncritical shock",
    "maximal_noncritical_shock_magnitude"   => "Max noncritical shock",
    "characteristic_return_time"            => "Char. return time",
    "reactivity"                            => "Reactivity",
    "maximal_amplification"                 => "Max amplification",
    "maximal_amplification_time"            => "Max amplification time",
    "intermingledness"                      => "Intermingledness",
)

# ─────────────────────────────────────────────────────────────────────────────
# Build the multi-parameter curve
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_co2_pcurve(p1x, p2x, t_values) → Vector{Vector{Pair{Int,Float64}}}

Return a parameter curve for `global_continuation`.  Each element is a list of
`(index => value)` pairs covering every parameter that differs between `p1x`
and `p2x`.  At step t the parameter vector is:

    params(t) = p1x + t * (p2x - p1x)
"""
function make_co2_pcurve(p1x, p2x, t_values)
    # Which parameter indices actually change?
    changing = findall(i -> p1x[i] != p2x[i], eachindex(p1x))

    pcurve = Vector{Vector{Pair{Int,Float64}}}(undef, length(t_values))
    for (k, t) in enumerate(t_values)
        p_t = p1x .+ t .* (p2x .- p1x)
        pcurve[k] = [i => p_t[i] for i in changing]
    end
    return pcurve
end

const T_VALUES   = collect(range(T_START, T_END; step = T_STEP))
const CO2_1X     = 280.0   # ppm, pre-industrial
const CO2_2X     = 560.0   # ppm, 2×CO₂
const CO2_VALUES = CO2_1X .+ T_VALUES .* (CO2_2X - CO2_1X)
const pcurve     = make_co2_pcurve(amoc_params_1xco2, amoc_params_2xco2, T_VALUES)

# ─────────────────────────────────────────────────────────────────────────────
# Dynamical system (base = 1xCO2; continuation updates params from there)
# ─────────────────────────────────────────────────────────────────────────────

ds, grid, proximity = amoc3box_1xco2()

sampler, = statespace_sampler(grid)

recurrence_kw = Dict(
    :Ttr                         => 0.0,
    :Δt                          => 1.0,
    :stop_at_Δt                  => true,
    :horizon_limit               => 1e2,
    :consecutive_lost_steps      => 10_000,
    :consecutive_recurrences     => 10_000,
    :consecutive_attractor_steps => 1_000,
    :consecutive_basin_steps     => 1_000,
)

proximity_options = (;
    Ttr                    = 0.0,
    Δt                     = 1.0,
    stop_at_Δt             = true,
    horizon_limit          = 1e2,
    consecutive_lost_steps = 10_000,
)

# ─────────────────────────────────────────────────────────────────────────────
# Identify AMOC-on attractor ID (highest overturning strength at t=0)
# ─────────────────────────────────────────────────────────────────────────────

function find_on_id(attractors_cont, params_base)
    for atts in attractors_cont
        isempty(atts) && continue
        id_q = [(id, amoc_strength(vec(mean(Matrix(att); dims = 1)), params_base))
                for (id, att) in atts if !isempty(att)]
        isempty(id_q) && continue
        return id_q[argmax(last.(id_q))][1]
    end
    return 1  # fallback
end

# ─────────────────────────────────────────────────────────────────────────────
# Run (or load cached) continuation + stability measures
# ─────────────────────────────────────────────────────────────────────────────

config = Dict(
    "T_start"              => T_START,
    "T_end"                => T_END,
    "T_step"               => T_STEP,
    "continuation_samples" => CONTINUATION_SAMPLES,
    "resilience_samples"   => RESILIENCE_SAMPLES,
    "finite_time"          => FINITE_TIME,
    "compute_stability"    => COMPUTE_STABILITY,
)

data, _ = produce_or_load(
    datadir("co2_continuation"), config;
    suffix = "jld2",
) do _

    mapper = AttractorsViaRecurrences(ds, grid; recurrence_kw...)
    ascm   = AttractorSeedContinueMatch(mapper)

    @info "Running CO2-trajectory global continuation..."
    fractions_cont, attractors_cont = global_continuation(
        ascm, pcurve, sampler;
        samples_per_parameter = CONTINUATION_SAMPLES,
    )

    # Identify AMOC-on attractor ID here so it is available inside produce_or_load
    on_id_inner = find_on_id(attractors_cont, copy(amoc_params_1xco2))

    nls_measures = if COMPUTE_STABILITY
        @info "Computing stability measures along continuation (attractor-centred sampling)..."
        stability_measures_attractor_centered(
            ds, attractors_cont, pcurve, on_id_inner;
            grid_delta               = 0.2,   # ±2 psu in model units
            grid_step                = 0.02,  # same resolution as the global grid
            ε                        = proximity.ε,
            weighting_distribution   = Attractors.EverywhereUniform(),
            finite_time              = FINITE_TIME,
            distance                 = StrictlyMinimumDistance(),
            samples_per_parameter    = RESILIENCE_SAMPLES,
            proximity_mapper_options = proximity_options,
        )
    else
        Dict{String, Any}()
    end

    @strdict fractions_cont attractors_cont nls_measures
end

fractions_cont  = data["fractions_cont"]
attractors_cont = data["attractors_cont"]
nls_measures    = data["nls_measures"]
n_steps         = length(T_VALUES)

on_id = find_on_id(attractors_cont, copy(amoc_params_1xco2))

# Last step at which the on-state attractor is still present
last_on = something(
    findlast(i -> haskey(attractors_cont[i], on_id) &&
                  !isempty(attractors_cont[i][on_id]), 1:n_steps),
    n_steps,
)

# ─────────────────────────────────────────────────────────────────────────────
# AMOC strength curve along t
# ─────────────────────────────────────────────────────────────────────────────

function amoc_q_along_pcurve(attractors_cont, pcurve, params_base)
    qs    = fill(NaN, length(attractors_cont))
    on_id = nothing

    for (i, (atts, updates)) in enumerate(zip(attractors_cont, pcurve))
        isempty(atts) && continue
        p = copy(params_base)
        for (idx, val) in updates
            p[idx] = val
        end

        id_q = [(id, amoc_strength(vec(mean(Matrix(att); dims = 1)), p))
                for (id, att) in atts if !isempty(att)]
        isempty(id_q) && continue

        if on_id === nothing
            on_id = id_q[argmax(last.(id_q))][1]
        end

        if any(id == on_id for (id, _) in id_q)
            qs[i] = first(q for (id, q) in id_q if id == on_id)
        else
            qs[i] = id_q[argmax(last.(id_q))][2]
        end
    end
    return qs
end

q_curve = amoc_q_along_pcurve(attractors_cont, pcurve, copy(amoc_params_1xco2))

# ─────────────────────────────────────────────────────────────────────────────
# Helper: extract one measure curve
# ─────────────────────────────────────────────────────────────────────────────

function measure_curve(nls_measures, measure_name, attractor_id, n_steps)
    out = fill(NaN, n_steps)
    !haskey(nls_measures, measure_name) && return out
    per_step = nls_measures[measure_name]
    for (i, d) in enumerate(per_step)
        i > n_steps && break
        haskey(d, attractor_id) && (out[i] = d[attractor_id])
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# Figure
# ─────────────────────────────────────────────────────────────────────────────

all_measures = COMPUTE_STABILITY ? sort(collect(keys(nls_measures))) : String[]
n_rows = length(all_measures) + 1

fig = Figure(size = (520, 180 * n_rows))

# ── Panel 1: AMOC overturning strength ──────────────────────────────────────
ax_q = Axis(fig[1, 1];
    ylabel             = "AMOC strength q (Sv)",
    xlabel             = isempty(all_measures) ? "CO₂ concentration (ppm)" : "",
    xticklabelsvisible = isempty(all_measures),
)

vq = .!isnan.(q_curve)
lines!(ax_q, CO2_VALUES[vq], q_curve[vq]; color = :steelblue, linewidth = 2)
hlines!(ax_q, [0.0]; color = :black, linestyle = :dash, linewidth = 1)

# Mark 1xCO2 and 2xCO2 positions
vlines!(ax_q, [CO2_1X, CO2_2X]; color = :gray, linestyle = :dot, linewidth = 1)
text!(ax_q, 0.01, 0.95; text = "1×CO₂\n(280 ppm)", align = (:left, :top),
      space = :relative, fontsize = 12, color = :gray)
text!(ax_q, (CO2_2X - CO2_1X) / (CO2_VALUES[end] - CO2_1X) + 0.01, 0.95;
      text = "2×CO₂\n(560 ppm)", align = (:left, :top),
      space = :relative, fontsize = 12, color = :gray)
text!(ax_q, 0.97, 0.05; text = "(a)", align = (:right, :bottom),
      space = :relative, fontsize = 16)

all_axes = Axis[ax_q]

# ── Panels 2…n_rows: stability measures (AMOC-on, clipped at tipping) ───────
if COMPUTE_STABILITY
    for (k, mname) in enumerate(all_measures)
        row     = k + 1
        is_last = row == n_rows
        ax = Axis(fig[row, 1];
            xlabel             = is_last ? "CO₂ concentration (ppm)" : "",
            ylabel             = get(MEASURE_YLABELS, mname, mname),
            xticklabelsvisible = is_last,
        )
        push!(all_axes, ax)

        cv = measure_curve(nls_measures, mname, on_id, n_steps)
        valid = [!isnan(cv[i]) && i <= last_on for i in 1:n_steps]

        lines!(ax, CO2_VALUES[valid], cv[valid]; color = :steelblue, linewidth = 2)
        vlines!(ax, [CO2_1X, CO2_2X]; color = :gray, linestyle = :dot, linewidth = 1)

        label = "($(('a':'z')[row]))"
        text!(ax, 0.97, 0.05; text = label, align = (:right, :bottom),
              space = :relative, fontsize = 16)
    end
end

linkxaxes!(all_axes...)

fig_path = plotsdir("amoc3box_co2_continuation.png")
wsave(fig_path, fig)
@info "Figure saved to: $fig_path"
display(fig)

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────

println("\n=== CO2 Trajectory Continuation Summary ===")
println("t range: $T_START → $T_END, step $T_STEP")
println("  t=0 : 1xCO2 (pre-industrial)")
println("  t=1 : 2xCO2 (~2°C warmer)")
println("Continuation samples : $CONTINUATION_SAMPLES per step")
println("Resilience samples   : $(COMPUTE_STABILITY ? RESILIENCE_SAMPLES : "skipped")")
println()

if COMPUTE_STABILITY
    bs = measure_curve(nls_measures, "basin_stability", on_id, n_steps)
    valid_bs = findall(!isnan, bs)
    if !isempty(valid_bs)
        zero_idx = findfirst(v -> v < 0.01, bs[valid_bs])
        if !isnothing(zero_idx)
            t_tip = T_VALUES[valid_bs[zero_idx]]
            println("AMOC-on basin stability → 0 near t ≈ $t_tip")
            println("  (1xCO2=0, 2xCO2=1 → $(round((t_tip-T_START)/(1.0-T_START)*100; digits=1))% along 1x→2x trajectory)")
        else
            println("AMOC-on attractor persists across full parameter range.")
        end
    end
end
println("AMOC-on last present at step $(last_on) / $(n_steps) (t = $(T_VALUES[last_on]))")

# ─────────────────────────────────────────────────────────────────────────────
# Basin plots along the CO₂ parameter curve
# ─────────────────────────────────────────────────────────────────────────────

if PLOT_BASINS
    basins_outdir = plotsdir("co2_continuation_basins")
    mkpath(basins_outdir)
    @info "Computing and plotting basins of attraction along CO₂ pcurve..."

    last_nonempty = findlast(atts -> !isempty(atts), attractors_cont)

    for (i, (updates, atts)) in enumerate(zip(pcurve, attractors_cont))
        isempty(atts) && continue

        # Apply parameter updates to ds for this pcurve step
        for (idx, val) in updates
            set_parameter!(ds, idx, val)
        end

        # Proximity mapper seeded with the attractors already found by continuation
        prox_mapper = AttractorsViaProximity(ds, atts;
            ε                      = proximity.ε,
            proximity_options...,
        )

        # Compute basins over the full state-space grid
        basins, attractors_step = basins_of_attraction(prox_mapper, grid;
            show_progress = true,
        )

        # Plot
        t_val   = T_VALUES[i]
        co2_val = CO2_VALUES[i]
        fig_b = heatmap_basins_attractors(grid, basins, attractors_step)
        ax_b  = fig_b.content[1]   # first Axis in the Figure
        ax_b.title  = "CO₂ = $(round(co2_val; digits=0)) ppm"
        ax_b.xlabel = "S_N (salinity anomaly ×100)"
        ax_b.ylabel = "S_T (salinity anomaly ×100)"

        step_str = lpad(string(i), 3, '0')
        t_str    = replace(string(round(t_val; digits=1)), "." => "p")
        fname    = joinpath(basins_outdir, "basins_$(step_str)_t$(t_str).png")
        wsave(fname, fig_b)
        @info "  Saved: $(basename(fname))"
    end

    # Reset ds parameters to 1×CO₂ baseline
    for (idx, val) in enumerate(amoc_params_1xco2)
        set_parameter!(ds, idx, val)
    end

    @info "Basin plots saved to: $basins_outdir"
end

# ─────────────────────────────────────────────────────────────────────────────
# Export resilience metrics to CSV for synthesis figure
# ─────────────────────────────────────────────────────────────────────────────
using CSV, DataFrames

rows = DataFrame(
    co2_ppm             = Float64[],
    t_param             = Float64[],
    measure             = String[],
    value               = Float64[],
    attractor           = String[],
)

# AMOC strength curve (convert m³/s → Sv: 1 Sv = 1e6 m³/s)
for (i, co2) in enumerate(CO2_VALUES)
    if !isnan(q_curve[i])
        push!(rows, (co2, T_VALUES[i], "amoc_strength_sv", q_curve[i] / 1e6, "on"))
    end
end

# Resilience measures — only the four used in the paper synthesis figure
const CSV_MEASURES = [
    "characteristic_return_time",
    "mean_convergence_time",
    "basin_stability",
    "minimal_critical_shock_magnitude",
]

if COMPUTE_STABILITY
    for mname in CSV_MEASURES
        cv = measure_curve(nls_measures, mname, on_id, n_steps)
        for (i, co2) in enumerate(CO2_VALUES)
            if !isnan(cv[i]) && i <= last_on
                push!(rows, (co2, T_VALUES[i], mname, cv[i], "on"))
            end
        end
    end
end

csv_path = datadir("paper", "resilience_vs_co2_boxmodel.csv")
mkpath(datadir("paper"))
CSV.write(csv_path, rows)
@info "Resilience CSV saved to: $csv_path"
