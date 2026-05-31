"""
    amoc3box_export_paper_data.jl

Export basin-of-attraction data, attractor positions, and representative
trajectories for the 3-box AMOC model to CSV files consumed by plotting_paper.py.

Loops over every t step that was computed in the CO₂ continuation and exports
one set of CSVs per step, tagged by the rounded CO₂ value:

    data/paper/basin_{co2}ppm.csv         — 60×60 basin label grid
    data/paper/trajectories_{co2}ppm.csv  — traj_id, time, S_N, S_T, q
    data/paper/attractors_{co2}ppm.csv    — state, S_N, S_T, q

Attractor locations are taken directly from the continuation cache (no fresh
attractor search per step).  Basin grids and trajectories are computed fresh.

Must be run after amoc3box_co2_continuation.jl has produced the JLD2 cache.

Run from the project root:
    julia --project scripts/amoc3box_export_paper_data.jl
"""

using DrWatson
@quickactivate "AMOCResilience"

include(srcdir("used_dynamical_systems.jl"))

using DynamicalSystems
using Attractors
using Statistics
using CSV
using DataFrames

# ─────────────────────────────────────────────────────────────────────────────
# Continuation parameters — must match amoc3box_co2_continuation.jl exactly
# so that produce_or_load finds the right JLD2 cache file.
# ─────────────────────────────────────────────────────────────────────────────

const T_START              = 0.0
const T_END                = 2.4
const T_STEP               = 0.05
const CONTINUATION_SAMPLES = 100
const RESILIENCE_SAMPLES   = 10_000
const FINITE_TIME          = 1000.0
const COMPUTE_STABILITY    = true

const CO2_1X = 280.0
const CO2_2X = 560.0

const T_VALUES   = collect(range(T_START, T_END; step = T_STEP))
const CO2_VALUES = CO2_1X .+ T_VALUES .* (CO2_2X - CO2_1X)

# ─────────────────────────────────────────────────────────────────────────────
# Settings
# ─────────────────────────────────────────────────────────────────────────────

const T_MAX      = 3000.0   # integration length for trajectories
const DT_TRAJ    = 1.0      # output time step
const GRID_RES   = 60       # basin grid resolution (60×60)
const BASIN_DELTA = 0.2     # half-width of attractor-centred basin grid (±2 psu)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: AMOC overturning strength in Sv
# ─────────────────────────────────────────────────────────────────────────────

sv(q_m3s) = q_m3s / 1e6

# ─────────────────────────────────────────────────────────────────────────────
# Helper: integrate ds from u0 for T_MAX years
# ─────────────────────────────────────────────────────────────────────────────

function integrate_trajectory(ds, u0, params, T_max, dt)
    reinit!(ds, u0)
    traj, tvec = DynamicalSystems.trajectory(ds, T_max; Δt = dt)
    n = length(tvec)
    S_N_arr = [traj[i][1] for i in 1:n]
    S_T_arr = [traj[i][2] for i in 1:n]
    q_arr   = [sv(amoc_strength([S_N_arr[i], S_T_arr[i]], params)) for i in 1:n]
    return tvec, S_N_arr, S_T_arr, q_arr
end

# ─────────────────────────────────────────────────────────────────────────────
# Helper: identify on-attractor ID (highest q)
# ─────────────────────────────────────────────────────────────────────────────

function find_on_id(attractors_dict, params)
    id_q = [(id, sv(amoc_strength(vec(mean(Matrix(att); dims = 1)), params)))
            for (id, att) in attractors_dict if !isempty(att)]
    isempty(id_q) && return nothing
    return id_q[argmax(last.(id_q))][1]
end

# ─────────────────────────────────────────────────────────────────────────────
# Process one t step
# ─────────────────────────────────────────────────────────────────────────────

function process_step(i, t_val, co2_val, attractors_from_cont)
    tag = "$(round(Int, co2_val))ppm"
    @info "Processing t=$(round(t_val; digits=2))  CO₂=$(round(co2_val; digits=0)) ppm  →  $tag"

    params = amoc_params_1xco2 .+ t_val .* (amoc_params_2xco2 .- amoc_params_1xco2)
    ds, grid, prox_kw = amoc3box_at_t(t_val)

    # ── Identify on- and off-attractors ──────────────────────────────────────
    on_id  = find_on_id(attractors_from_cont, params)
    if isnothing(on_id)
        @warn "  No valid attractor found at step $i (t=$t_val) — skipping."
        return
    end

    id_q_pairs = [
        (id, sv(amoc_strength(vec(mean(Matrix(att); dims = 1)), params)))
        for (id, att) in attractors_from_cont if !isempty(att)
    ]
    sort!(id_q_pairs; by = last, rev = true)
    off_id = length(id_q_pairs) >= 2 ? id_q_pairs[2][1] : on_id

    att_on_state  = vec(mean(Matrix(attractors_from_cont[on_id]);  dims = 1))
    att_off_state = (on_id != off_id) ?
                    vec(mean(Matrix(attractors_from_cont[off_id]); dims = 1)) :
                    att_on_state .+ [0.5, -0.5]  # fallback: monostable step

    q_on  = sv(amoc_strength(att_on_state,  params))
    q_off = sv(amoc_strength(att_off_state, params))

    # ── Export attractor positions ────────────────────────────────────────────
    attr_df = DataFrame(
        state = ["on",            "off"           ],
        S_N   = [att_on_state[1], att_off_state[1]],
        S_T   = [att_on_state[2], att_off_state[2]],
        q     = [q_on,            q_off           ],
    )
    attr_path = datadir("paper", "attractors_$(tag).csv")
    CSV.write(attr_path, attr_df)
    @info "  Attractor CSV saved: $attr_path"

    # ── Basin grid centred on on-attractor ────────────────────────────────────
    sn_c, st_c = att_on_state[1], att_on_state[2]
    xg = range(sn_c - BASIN_DELTA, sn_c + BASIN_DELTA; length = GRID_RES)
    yg = range(st_c - BASIN_DELTA, st_c + BASIN_DELTA; length = GRID_RES)

    prox_mapper = AttractorsViaProximity(ds, attractors_from_cont;
        ε = prox_kw.ε,
        Ttr                    = prox_kw.Ttr,
        Δt                     = prox_kw.Δt,
        stop_at_Δt             = prox_kw.stop_at_Δt,
        horizon_limit          = prox_kw.horizon_limit,
        consecutive_lost_steps = prox_kw.consecutive_lost_steps,
    )
    basins, _ = basins_of_attraction(prox_mapper, (xg, yg); show_progress = false)

    basin_rows = DataFrame(S_N = Float64[], S_T = Float64[], basin_label = Int[])
    for (ix, sn) in enumerate(xg), (iy, st) in enumerate(yg)
        push!(basin_rows, (sn, st, basins[ix, iy]))
    end
    basin_path = datadir("paper", "basin_$(tag).csv")
    CSV.write(basin_path, basin_rows)
    @info "  Basin CSV saved: $basin_path (centred on on-attractor, ±$(BASIN_DELTA) model units)"

    # ── Trajectories from four corners of the attractor-centred grid ─────────
    r = BASIN_DELTA * 0.8   # slightly inside the grid boundary
    ics = [
        SVector(sn_c - r, st_c + r),   # upper-left
        SVector(sn_c + r, st_c + r),   # upper-right
        SVector(sn_c - r, st_c - r),   # lower-left
        SVector(sn_c + r, st_c - r),   # lower-right
    ]

    traj_rows = DataFrame(
        traj_id = Int[],
        time    = Float64[],
        S_N     = Float64[],
        S_T     = Float64[],
        q       = Float64[],
    )

    for (tid, ic) in enumerate(ics)
        @info "  Integrating trajectory $tid/4 from $ic ..."
        tvec, SN_arr, ST_arr, q_arr = integrate_trajectory(ds, ic, params, T_MAX, DT_TRAJ)
        for k in 1:length(tvec)
            push!(traj_rows, (tid, tvec[k], SN_arr[k], ST_arr[k], q_arr[k]))
        end
    end

    traj_path = datadir("paper", "trajectories_$(tag).csv")
    CSV.write(traj_path, traj_rows)
    @info "  Trajectories CSV saved: $traj_path"

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

function main()
    mkpath(datadir("paper"))

    # ── Load continuation results from cache ─────────────────────────────────
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
        error(
            "Continuation cache not found. " *
            "Run amoc3box_co2_continuation.jl first."
        )
    end

    attractors_cont = data["attractors_cont"]
    n_steps = length(T_VALUES)

    # Identify on-attractor ID using the 1×CO₂ parameters
    on_id_global = let
        first_nonempty = findfirst(!isempty, attractors_cont)
        isnothing(first_nonempty) && error("No attractors found in continuation.")
        find_on_id(attractors_cont[first_nonempty], copy(amoc_params_1xco2))
    end
    @info "On-attractor ID: $on_id_global"

    # ── Export every step where the on-attractor is present ──────────────────
    n_exported = 0
    for i in 1:n_steps
        atts = attractors_cont[i]
        isempty(atts)                          && continue
        !haskey(atts, on_id_global)            && continue
        isempty(atts[on_id_global])            && continue

        process_step(i, T_VALUES[i], CO2_VALUES[i], atts)
        n_exported += 1
    end

    @info "Done. Exported $n_exported steps to $(datadir("paper"))"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
