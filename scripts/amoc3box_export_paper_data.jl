"""
    amoc3box_export_paper_data.jl

Export basin-of-attraction data, attractor positions, and representative
trajectories for the 3-box AMOC model to CSV files consumed by plotting_paper.py.

Outputs (all in data/paper/):
    basin_{scenario}.csv         — S_N, S_T, basin_label columns (60×60 grid)
    trajectories_{scenario}.csv  — traj_id, time, S_N, S_T, q
    attractors_{scenario}.csv    — state, S_N, S_T, q

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
# Settings
# ─────────────────────────────────────────────────────────────────────────────

const T_MAX      = 3000.0   # integration length (model years) for IC3, IC4
const T_MAX_LONG = 6000.0   # integration length for IC1, IC2
const DT_TRAJ    = 1.0      # output time step
const GRID_RES   = 60       # basin grid resolution (60×60)
const BASIN_DELTA = 0.2     # half-width of attractor-centred basin grid (±2 psu)

# ─────────────────────────────────────────────────────────────────────────────
# Helper: AMOC overturning strength in Sv
#   q comes out in m³/s from amoc_strength(); 1 Sv = 1e6 m³/s
# ─────────────────────────────────────────────────────────────────────────────

sv(q_m3s) = q_m3s / 1e6

# ─────────────────────────────────────────────────────────────────────────────
# Helper: integrate ds from u0 for T_MAX years, recording (t, S_N, S_T, q)
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
# Process one scenario
# ─────────────────────────────────────────────────────────────────────────────

function process_scenario(scenario::String)
    @info "Processing scenario: $scenario"

    # Build dynamical system
    ds, grid_full, proximity = if scenario == "1xco2"
        amoc3box_1xco2()
    elseif scenario == "2xco2"
        amoc3box_2xco2()
    else
        amoc3box_896ppm()
    end
    params = if scenario == "1xco2"
        copy(amoc_params_1xco2)
    elseif scenario == "2xco2"
        copy(amoc_params_2xco2)
    else
        copy(amoc_params_896ppm)
    end

    # ── Step 1: find attractors on the full grid via recurrences ────────────────
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

    xg_full = range(grid_full[1][1], grid_full[1][end]; length = GRID_RES)
    yg_full = range(grid_full[2][1], grid_full[2][end]; length = GRID_RES)

    mapper_full = AttractorsViaRecurrences(ds, (xg_full, yg_full); recurrence_kw...)
    _, attractors = basins_of_attraction(mapper_full, (xg_full, yg_full); show_progress = true)

    # ── Identify on-attractor (highest q) and off-attractor ──────────────────
    id_q_pairs = [
        (id, sv(amoc_strength(vec(mean(Matrix(att); dims = 1)), params)))
        for (id, att) in attractors if !isempty(att)
    ]
    sort!(id_q_pairs; by = last, rev = true)

    on_id  = id_q_pairs[1][1]
    off_id = length(id_q_pairs) >= 2 ? id_q_pairs[2][1] : id_q_pairs[1][1]

    # ── Step 2: compute basin on ±BASIN_DELTA grid centred on on-attractor ───
    att_on_centroid = vec(mean(Matrix(attractors[on_id]); dims = 1))
    sn_c, st_c = att_on_centroid[1], att_on_centroid[2]

    xg = range(sn_c - BASIN_DELTA, sn_c + BASIN_DELTA; length = GRID_RES)
    yg = range(st_c - BASIN_DELTA, st_c + BASIN_DELTA; length = GRID_RES)

    proximity_kw = (; Ttr = 0.0, Δt = 1.0, stop_at_Δt = true,
                      horizon_limit = 1e2, consecutive_lost_steps = 10_000)
    prox_mapper = AttractorsViaProximity(ds, attractors;
        ε = proximity.ε, proximity_kw...)
    basins, _ = basins_of_attraction(prox_mapper, (xg, yg); show_progress = true)

    # ── Export basin CSV ─────────────────────────────────────────────────────
    basin_rows = DataFrame(S_N = Float64[], S_T = Float64[], basin_label = Int[])
    for (ix, sn) in enumerate(xg), (iy, st) in enumerate(yg)
        push!(basin_rows, (sn, st, basins[ix, iy]))
    end
    basin_path = datadir("paper", "basin_$(scenario).csv")
    CSV.write(basin_path, basin_rows)
    @info "  Basin CSV saved: $basin_path (centred on on-attractor, ±$(BASIN_DELTA) model units = ±$(Int(BASIN_DELTA*10)) psu)"

    att_on_state  = vec(mean(Matrix(attractors[on_id]);  dims = 1))
    att_off_state = (on_id != off_id) ? vec(mean(Matrix(attractors[off_id]); dims = 1)) :
                    att_on_state .+ [0.5, -0.5]  # fallback if monostable

    q_on  = sv(amoc_strength(att_on_state,  params))
    q_off = sv(amoc_strength(att_off_state, params))

    # ── Export attractor positions ────────────────────────────────────────────
    attr_df = DataFrame(
        state = ["on",         "off"        ],
        S_N   = [att_on_state[1],  att_off_state[1] ],
        S_T   = [att_on_state[2],  att_off_state[2] ],
        q     = [q_on,             q_off            ],
    )
    attr_path = datadir("paper", "attractors_$(scenario).csv")
    CSV.write(attr_path, attr_df)
    @info "  Attractor CSV saved: $attr_path"

    # ── Select 4 initial conditions (scenario-specific) ──────────────────────
    if scenario == "1xco2"
        ic1 = SVector(-0.2,  0.15)
        ic2 = SVector(-0.15, -0.2)   # S_T perturbation sign flipped
        ic3 = SVector( 0.1,  0.1 )
        ic4 = SVector( 0.2,  0.0 )
    else  # 896ppm — all ICs within the ±0.2 perturbation box around on-attractor (0.288, 0.804)
        ic1 = SVector( 0.45,  0.95)  # upper-right
        ic2 = SVector( 0.15,  0.95)  # upper-left
        ic3 = SVector( 0.45,  0.65)  # lower-right
        ic4 = SVector( 0.15,  0.65)  # lower-left
    end

    ics    = [ic1,        ic2,        ic3,    ic4   ]
    t_maxs = [T_MAX_LONG, T_MAX_LONG, T_MAX,  T_MAX ]

    # ── Integrate trajectories ────────────────────────────────────────────────
    traj_rows = DataFrame(
        traj_id = Int[],
        time    = Float64[],
        S_N     = Float64[],
        S_T     = Float64[],
        q       = Float64[],
    )

    for (tid, (ic, t_max)) in enumerate(zip(ics, t_maxs))
        @info "  Integrating trajectory $tid / 4 from $ic for $t_max years..."
        tvec, SN_arr, ST_arr, q_arr = integrate_trajectory(ds, ic, params, t_max, DT_TRAJ)
        for k in 1:length(tvec)
            push!(traj_rows, (tid, tvec[k], SN_arr[k], ST_arr[k], q_arr[k]))
        end
    end

    traj_path = datadir("paper", "trajectories_$(scenario).csv")
    CSV.write(traj_path, traj_rows)
    @info "  Trajectories CSV saved: $traj_path"

    return nothing
end

# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

function main()
    mkpath(datadir("paper"))

    for scenario in ("1xco2", "896ppm")
        process_scenario(scenario)
    end

    @info "All paper data exported to $(datadir("paper"))"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
