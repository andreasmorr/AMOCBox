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
const RESILIENCE_SAMPLES   = 1000   # samples per parameter step for resilience

# Finite-time horizon for finite-time basin stability
const FINITE_TIME = 1000.0

# Resilience measures to extract for plotting
const CHOSEN_MEASURES = [
    "minimal_critical_shock_magnitude",
    "maximal_noncritical_shock_magnitude",
    "basin_stability",
    "median_convergence_time",
    "finite_time_basin_stability",
]

const MEASURE_YLABELS = Dict(
    "minimal_critical_shock_magnitude"    => "Min. critical shock",
    "maximal_noncritical_shock_magnitude" => "Max. non-critical shock",
    "basin_stability"                     => "Basin stability",
    "median_convergence_time"             => "Med. convergence time (yr)",
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

    # Recurrences mapper keywords
    recurrence_kw = Dict(
        :ε                       => ε,
        :Ttr                     => proximity.Ttr,
        :Δt                      => proximity.Δt,
        :stop_at_Δt              => proximity.stop_at_Δt,
        :horizon_limit           => proximity.horizon_limit,
        :consecutive_lost_steps  => proximity.consecutive_lost_steps,
    )

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
        proximity_mapper_options   = proximity,
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
    qs = fill(NaN, length(H_values))
    for (i, (h, atts)) in enumerate(zip(H_values, attractors_cont))
        p = copy(params_base)
        p[1] = h
        # Find the on-attractor (highest q) among all attractors at this H
        best_q = -Inf
        for (id, att) in atts
            isempty(att) && continue
            mean_state = vec(mean(Matrix(att); dims = 1))
            q = amoc_strength(mean_state, p)
            best_q = max(best_q, q)
        end
        qs[i] = isinf(best_q) ? NaN : best_q
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
fig = Figure(size = (800, 300 * n_rows))

for (row, mname) in enumerate(CHOSEN_MEASURES)
    ax = Axis(fig[row, 1];
        xlabel    = row == length(CHOSEN_MEASURES) ? "Hosing H (Sv)" : "",
        ylabel    = get(MEASURE_YLABELS, mname, mname),
        xticklabelsvisible = row == length(CHOSEN_MEASURES),
    )

    c1 = measure_curve(nls_1x, mname, ID_ON, n_steps)
    c2 = measure_curve(nls_2x, mname, ID_ON, n_steps)

    v1 = .!isnan.(c1);  v2 = .!isnan.(c2)

    lines!(ax, H_values[v1], c1[v1]; color = col_1x, linewidth = 2, label = "1×CO₂")
    lines!(ax, H_values[v2], c2[v2]; color = col_2x, linewidth = 2, label = "2×CO₂")

    row == 1 && axislegend(ax; position = :rt)

    # Letter label
    text!(ax, 0.97, 0.05;
        text       = "($(('a':'z')[row]))",
        align      = (:right, :bottom),
        space      = :relative,
        fontsize   = 16,
    )
end

# Bottom panel: AMOC overturning strength
ax_q = Axis(fig[n_rows, 1];
    xlabel = "Hosing H (Sv)",
    ylabel = "AMOC strength q (m³/s)",
)
vq1 = .!isnan.(q_1x);  vq2 = .!isnan.(q_2x)
lines!(ax_q, H_values[vq1], q_1x[vq1]; color = col_1x, linewidth = 2, label = "1×CO₂")
lines!(ax_q, H_values[vq2], q_2x[vq2]; color = col_2x, linewidth = 2, label = "2×CO₂")
hlines!(ax_q, [0.0]; color = :black, linestyle = :dash, linewidth = 1)
axislegend(ax_q; position = :rt)
text!(ax_q, 0.97, 0.05;
    text = "($(('a':'z')[n_rows]))", align = (:right, :bottom),
    space = :relative, fontsize = 16,
)

resize!(fig, 700, 280 * n_rows)
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
