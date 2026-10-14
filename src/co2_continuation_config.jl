"""
    co2_continuation_config.jl

Shared definition of the CO₂ continuation used by

    scripts/amoc3box_co2_continuation.jl   (runs the continuation, writes the JLD2 cache)
    scripts/amoc3box_export_paper_data.jl  (reads the JLD2 cache)

Both scripts include this file so that the parameter path and the
`produce_or_load` configuration (and hence the cache file name) cannot drift
apart.

Parameter path:  params(t) = params_1x + t * (params_2x - params_1x)
CO₂ label:       CO₂(t)    = 280 ppm + 280 ppm * t
"""

const T_START = 0.0    # 1xCO2 (280 ppm)
const T_END   = 2.4    # extrapolation beyond 2xCO2 (952 ppm)
const T_STEP  = 0.05   # regular spacing (14 ppm)

# Extra continuation points inserted between the regular steps.
# t = 2.225 (903 ppm) resolves the last AMOC-on step between 896 ppm (t = 2.20)
# and 910 ppm (t = 2.25), where the on-state is no longer found.
const T_REFINE = [2.225]

const CONTINUATION_SAMPLES = 100
const RESILIENCE_SAMPLES   = 10_000
const FINITE_TIME          = 1000.0
const COMPUTE_STABILITY    = true   # set false to skip stability measures (faster)

# Base seed for all random initial-condition sampling. The global-continuation
# sampler uses SEED; the attractor-centred resilience sampler at continuation
# step i uses SEED + i. With a fixed seed the whole run is reproducible
# (all mappers used here are single-threaded).
const SEED = 1234

# Maximum number of Δt = 1 yr steps an initial condition may take to come within
# ε of an attractor before AttractorsViaProximity labels it unassigned (-1).
# Near the end of the on-branch the on-state decays very slowly (e-folding time
# ≈ 3000 yr at 903 ppm) and convergence from the ±2 psu box takes up to
# ≈ 18,000 yr; the Attractors.jl default of 10,000 left ~25% of the samples at
# 903 ppm unassigned. Used for the resilience measures and the exported basins.
const CONSECUTIVE_LOST_STEPS = 50_000

const CO2_1X = 280.0   # ppm, pre-industrial
const CO2_2X = 560.0   # ppm, 2×CO₂

const T_VALUES   = sort(unique(vcat(collect(range(T_START, T_END; step = T_STEP)), T_REFINE)))
const CO2_VALUES = CO2_1X .+ T_VALUES .* (CO2_2X - CO2_1X)

"""
    continuation_config() → Dict

Configuration passed to DrWatson's `produce_or_load`. Every entry is a scalar
or string so that DrWatson includes it in the cache file name; in particular
`T_refine` makes the cache name change whenever the refinement points change.
"""
continuation_config() = Dict(
    "T_start"              => T_START,
    "T_end"                => T_END,
    "T_step"               => T_STEP,
    "T_refine"             => join(string.(T_REFINE), ","),
    "continuation_samples" => CONTINUATION_SAMPLES,
    "resilience_samples"   => RESILIENCE_SAMPLES,
    "finite_time"          => FINITE_TIME,
    "compute_stability"    => COMPUTE_STABILITY,
    "seed"                 => SEED,
    "consecutive_lost_steps" => CONSECUTIVE_LOST_STEPS,
)
