"""
    amoc3box_hosing_scan.jl

Hosing parameter (H) continuation for the 3-box AMOC model at 1xCO2 and 2xCO2.

For each CO2 scenario the script:
1. Sets up the dynamical system using `amoc3box_1xco2()` / `amoc3box_2xco2()`.
2. Runs global attractor continuation over H = 0 → 0.55 Sv using
   `AttractorSeedContinueMatch` + `global_continuation`.
3. Computes stability / resilience measures along the continuation using
   `stability_measures_along_continuation`.
4. Caches results via DrWatson's `produce_or_load`.
5. Creates a multi-panel comparison figure and saves it to `plotsdir()`.

Run from the project root:
    julia --project scripts/amoc3box_hosing_scan.jl
"""

using DrWatson
@quickactivate "AMOCResilience"

include(srcdir("used_dynamical_systems.jl"))

using DynamicalSystems
using Attractors
using CairoMakie
using Statistics

# ─────────────────────────────────────────────────────────────────────────────
# Tunable parameters
# ─────────────────────────────────────────────────────────────────────────────

const H_START  = 0.0
const H_END    = 0.55
const H_STEP   = 0.01

# Samples for basin / continuation map
const CONTINUATION_SAMPLES = 100    # samples per parameter step for continuation
const RESILIENCE_SAMPLES   = 10_000   # samples per parameter step for resilience

# Finite-time horizon for finite-time basin stability
const FINITE_TIME = 1000.0

# Resilience measures to extract for plotting
const CHOSEN_MEASURES = [
    "minimal_critical_shock_magnitude",
    "basin_stability",
    "median_convergence_time",
    "median_convergence_pace",
    "finite_time_basin_stability",
]

const MEASURE_YLABELS = Dict(
    "minimal_critical_shock_magnitude"    => "Min. critical shock",
    "basin_stability"                     => "Basin stability",
    "median_convergence_time"             => "Med. convergence time (yr)",
    "median_convergence_pace"             => "Med. convergence pace",
    "finite_time_basin_stability"         => "Finite-time basin stab.",
)

# Attractor IDs (as resolved by the matcher; 1 = AMOC-on, 2 = AMOC-off)
# You may need to swap these if the matcher assigns IDs differently.
const ID_ON  = 1
const ID_OFF = 2

# ─────────────────────────────────────────────────────────────────────────────
# Helper: run continuation + resilience for one CO2 scenario
# ─────────────────────────────────────────────────────────────────────────────

"""
    run_continuation(ds_fn, scenario_label; force=false) → Dict

Run global attractor continuation over the hosing parameter H and compute
stability measures along the continuation for the dynamical system returned
by the zero-argument function `ds_fn`.

Results are cached via DrWatson's `produce_or_load`.
"""
function run_continuation(ds_fn::Function, scenario_label::String;
                          force::Bool = false)
    config = Dict(
        "scenario"             => scenario_label,
        "H_start"              => H_START,
        "H_end"                => H_END,
        "H_step"               => H_STEP,
        "continuation_samples" => CONTINUATION_SAMPLES,
        "resilience_samples"   => RESILIENCE_SAMPLES,
        "finite_time"          => FINITE_TIME,
    )

    data, _ = produce_or_load(
        datadir("hosing_continuation"), config;
        suffix = "jld2", force,
    ) do cfg
        _continuation_inner(ds_fn, cfg)
    end

    return data
end

function _continuation_inner(ds_fn::Function, cfg::Dict)
    pstart = cfg["H_start"]
    pend   = cfg["H_end"]
    pstep  = cfg["H_step"]

    ds, grid, proximity = ds_fn()
    ε = proximity.ε

    # Parameter curve: each element is a list of parameter changes (index => value)
    # H is parameter index 1 in the parameter vector
    pcurve = [[1 => H] for H in range(pstart, pend; step = pstep)]

    # Sampler over state space
    sampler, = statespace_sampler(grid)

    # Recurrences mapper keywords (ε removed: current API uses grid cell size for proximity)
    recurrence_kw = Dict(
        :Ttr                     => proximity.Ttr,
        :Δt                      => proximity.Δt,
        :stop_at_Δt              => proximity.stop_at_Δt,
        :horizon_limit           => proximity.horizon_limit,
        :consecutive_lost_steps  => proximity.consecutive_lost_steps,
        :consecutive_recurrences => 10000,
        :consecutive_attractor_steps => 1000,   
        :consecutive_basin_steps => 1000,
    )
    # proximity_mapper_options for AttractorsViaProximity (used in stability measures)
    # must not include ε since it is passed as a separate keyword
    proximity_options = (; Ttr = proximity.Ttr, Δt = proximity.Δt,
                           stop_at_Δt = proximity.stop_at_Δt,
                           horizon_limit = proximity.horizon_limit,
                           consecutive_lost_steps = proximity.consecutive_lost_steps)

    mapper = AttractorsViaRecurrences(ds, grid; recurrence_kw...)
    ascm   = AttractorSeedContinueMatch(mapper)

    # Global continuation: track attractors as H increases
    fractions_cont, attractors_cont = global_continuation(
        ascm, pcurve, sampler;
        samples_per_parameter = cfg["continuation_samples"],
    )

    # Resilience / stability measures along the continuation
    nls_measures = stability_measures_along_continuation(
        ds, attractors_cont, pcurve, sampler;
        ε                          = ε,
        weighting_distribution     = Attractors.EverywhereUniform(),
        finite_time                = cfg["finite_time"],
        distance                   = StrictlyMinimumDistance(),
        samples_per_parameter      = cfg["resilience_samples"],
        proximity_mapper_options   = proximity_options,
    )

    return @strdict fractions_cont attractors_cont nls_measures pcurve
end

# ─────────────────────────────────────────────────────────────────────────────
# Run both scenarios
# ─────────────────────────────────────────────────────────────────────────────

@info "Running 1xCO2 continuation..."
data_1x = run_continuation(amoc3box_1xco2, "1xco2")

@info "Running 2xCO2 continuation..."
data_2x = run_continuation(amoc3box_2xco2, "2xco2")

# ─────────────────────────────────────────────────────────────────────────────
# Extract results
# ─────────────────────────────────────────────────────────────────────────────

H_values = [p[1][2] for p in data_1x["pcurve"]]   # extract H from pcurve
n_steps  = length(H_values)

nls_1x = data_1x["nls_measures"]
nls_2x = data_2x["nls_measures"]
attractors_cont_1x = data_1x["attractors_cont"]
attractors_cont_2x = data_2x["attractors_cont"]

# AMOC strength q at each H, for the on-attractor
function amoc_q_curve(attractors_cont, params_base, H_values)
    qs    = fill(NaN, length(H_values))
    on_id = nothing   # determined at first non-empty step: attractor with highest q

    for (i, (h, atts)) in enumerate(zip(H_values, attractors_cont))
        isempty(atts) && continue
        p = copy(params_base)
        p[1] = h

        # q for each non-empty attractor at this step
        id_q = [(id, amoc_strength(vec(mean(Matrix(att); dims=1)), p))
                for (id, att) in atts if !isempty(att)]
        isempty(id_q) && continue

        if on_id === nothing
            # First step: identify the on-state as the attractor with highest q
            on_id = id_q[argmax(last.(id_q))][1]
        end

        if any(id == on_id for (id, _) in id_q)
            # On-state still present — follow it
            qs[i] = first(q for (id, q) in id_q if id == on_id)
        else
            # On-state gone after tipping — plot whichever attractor remains
            qs[i] = id_q[argmax(last.(id_q))][2]
        end
    end
    return qs
end

q_1x = amoc_q_curve(attractors_cont_1x, copy(amoc_params_1xco2), H_values)
q_2x = amoc_q_curve(attractors_cont_2x, copy(amoc_params_2xco2), H_values)

# Helper: extract a resilience measure curve for a specific attractor ID
function measure_curve(nls_measures, measure_name::String, attractor_id::Int,
                       n_steps::Int)
    out = fill(NaN, n_steps)
    !haskey(nls_measures, measure_name) && return out
    per_step = nls_measures[measure_name]   # Vector of Dicts
    for (i, d) in enumerate(per_step)
        i > n_steps && break
        haskey(d, attractor_id) && (out[i] = d[attractor_id])
    end
    return out
end

# ─────────────────────────────────────────────────────────────────────────────
# Figure: multi-panel resilience comparison
# ─────────────────────────────────────────────────────────────────────────────

col_1x = :steelblue
col_2x = :firebrick

n_rows = length(CHOSEN_MEASURES) + 1   # +1 for AMOC strength panel

# Identify the on-state attractor ID at the first step (highest q = on-state)
function find_on_id(attractors_cont, params_base, H_values)
    for (h, atts) in zip(H_values, attractors_cont)
        isempty(atts) && continue
        p = copy(params_base); p[1] = h
        id_q = [(id, amoc_strength(vec(mean(Matrix(att); dims=1)), p))
                for (id, att) in atts if !isempty(att)]
        isempty(id_q) && continue
        return id_q[argmax(last.(id_q))][1]
    end
    return ID_ON  # fallback
end

on_id_1x = find_on_id(attractors_cont_1x, copy(amoc_params_1xco2), H_values)
on_id_2x = find_on_id(attractors_cont_2x, copy(amoc_params_2xco2), H_values)

# Last H step where AMOC-on genuinely exists — used to clip resilience curves
last_1x = something(
    findlast(i -> haskey(attractors_cont_1x[i], on_id_1x) && !isempty(attractors_cont_1x[i][on_id_1x]), 1:n_steps),
    n_steps)
last_2x = something(
    findlast(i -> haskey(attractors_cont_2x[i], on_id_2x) && !isempty(attractors_cont_2x[i][on_id_2x]), 1:n_steps),
    n_steps)

fig = Figure(size = (500, 180 * n_rows))

# Panel 1: AMOC overturning strength (on-state → off-state after tipping)
ax_q = Axis(fig[1, 1];
    ylabel = "AMOC strength q (Sv)",
    xticklabelsvisible = false,
)
vq1 = .!isnan.(q_1x);  vq2 = .!isnan.(q_2x)
lines!(ax_q, H_values[vq1], q_1x[vq1]; color = col_1x, linewidth = 2, label = "1×CO₂")
lines!(ax_q, H_values[vq2], q_2x[vq2]; color = col_2x, linewidth = 2, label = "2×CO₂")
hlines!(ax_q, [0.0]; color = :black, linestyle = :dash, linewidth = 1)
axislegend(ax_q; position = :rt)
text!(ax_q, 0.97, 0.05; text = "(a)", align = (:right, :bottom), space = :relative, fontsize = 16)

all_axes = Axis[ax_q]

# Panels 2..n_rows: resilience measures for AMOC-on, clipped at tipping point
for (k, mname) in enumerate(CHOSEN_MEASURES)
    row     = k + 1
    is_last = row == n_rows
    ax = Axis(fig[row, 1];
        xlabel             = is_last ? "Hosing H (Sv)" : "",
        ylabel             = get(MEASURE_YLABELS, mname, mname),
        xticklabelsvisible = is_last,
    )
    push!(all_axes, ax)

    c1 = measure_curve(nls_1x, mname, on_id_1x, n_steps)
    c2 = measure_curve(nls_2x, mname, on_id_2x, n_steps)

    v1 = [!isnan(c1[i]) && i <= last_1x for i in 1:n_steps]
    v2 = [!isnan(c2[i]) && i <= last_2x for i in 1:n_steps]

    lines!(ax, H_values[v1], c1[v1]; color = col_1x, linewidth = 2, label = "1×CO₂")
    lines!(ax, H_values[v2], c2[v2]; color = col_2x, linewidth = 2, label = "2×CO₂")

    text!(ax, 0.97, 0.05;
        text     = "($(('a':'z')[row]))",
        align    = (:right, :bottom),
        space    = :relative,
        fontsize = 16,
    )
end

linkxaxes!(all_axes...)

fig_path = plotsdir("amoc3box_hosing_resilience_comparison.png")
wsave(fig_path, fig)
@info "Figure saved to: $fig_path"
display(fig)

# ─────────────────────────────────────────────────────────────────────────────
# Summary printout
# ─────────────────────────────────────────────────────────────────────────────

println("\n=== AMOC 3-Box Hosing Scan Summary ===")
println("H range: $H_START — $H_END Sv, step $H_STEP")
println("Continuation samples : $CONTINUATION_SAMPLES per H step")
println("Resilience samples   : $RESILIENCE_SAMPLES per H step")
println()

for (label, nls) in [("1xCO2", nls_1x), ("2xCO2", nls_2x)]
    bs = measure_curve(nls, "basin_stability", ID_ON, n_steps)
    valid = findall(!isnan, bs)
    if !isempty(valid)
        zero_idx = findfirst(v -> v < 0.01, bs[valid])
        if !isnothing(zero_idx)
            h_tip = H_values[valid[zero_idx]]
            println("$label: AMOC-on basin stability → 0 near H ≈ $h_tip Sv")
        else
            println("$label: AMOC-on attractor persists across full H range scanned")
        end
    end
end
