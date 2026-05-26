"""
    used_dynamical_systems.jl

Defines the 3-box AMOC model (Wood 1993, Alkhayuon et al. 2019) for two CO2
scenarios: 1xCO2 (pre-industrial baseline) and 2xCO2 (~2°C warmer world).

Each scenario returns the ODE system, a state-space grid for attractor finding,
and a proximity threshold for attractor identification.

References:
  - Wood, R. (1993). Three-box model of AMOC. Internal UKMO report.
  - Alkhayuon, H. et al. (2019). Basin boundaries and bifurcations. Proc. R. Soc. A.
"""

using StaticArrays
using DynamicalSystems
using Attractors
using OrdinaryDiffEqRosenbrock

# ─────────────────────────────────────────────────────────────────────────────
# Parameter vectors
#
# Parameter order:
#   H      - freshwater hosing flux (Sv)
#   V_N    - volume North Atlantic box (m³)
#   V_T    - volume tropical box (m³)
#   V_S    - volume southern box (m³)
#   V_IP   - volume Indo-Pacific box (m³)
#   V_B    - volume deep box (m³)
#   S_S    - salinity southern box (psu)
#   S_IP   - salinity Indo-Pacific box (psu) [diagnostically determined]
#   S_B    - salinity deep box (psu)
#   F_N    - freshwater flux into North box (m³/s)  [before hosing]
#   F_T    - freshwater flux into tropical box (m³/s) [before hosing]
#   α      - thermal expansion coefficient (1/°C)
#   β      - haline contraction coefficient (1/psu)
#   S_0    - reference salinity (psu)
#   T_S    - southern surface temperature (°C) -- CHANGES between CO2 scenarios
#   T_0    - deep ocean temperature (°C)        -- CHANGES between CO2 scenarios
#   K_N    - mixing coefficient North (m³/s)
#   K_S    - mixing coefficient South (m³/s)
#   λ      - hydraulic constant (m³/s/Pa)
#   γ      - mixing fraction from southern/Indo-Pacific box
#   μ      - momentum advection coefficient
#   Y      - scaling factor (s/yr × unit conversions)
#   C      - total salinity content (conservation constraint)
# ─────────────────────────────────────────────────────────────────────────────

"""
    amoc_params_1xco2

Parameter vector for the 3-box AMOC model at 1xCO2.
"""
const amoc_params_1xco2 = [
    0.0,         # H:    hosing flux (Sv), start at zero
    0.3261e17,   # V_N:  North box volume (m³)
    0.7777e17,   # V_T:  tropical box volume (m³)
    0.8897e17,   # V_S:  southern box volume (m³)
    2.2020e17,   # V_IP: Indo-Pacific box volume (m³)
    8.6490e17,   # V_B:  deep box volume (m³)
    0.034427,    # S_S:  southern salinity (psu, dimensional)
    0.034668,    # S_IP: Indo-Pacific salinity (psu, initial guess; overwritten in ODE)
    0.034538,    # S_B:  deep salinity (psu)
    0.384e6,     # F_N:  freshwater flux North (m³/s)
    -0.723e6,    # F_T:  freshwater flux tropical (m³/s)
    0.12,        # α:    thermal expansion (1/°C)
    790.0,       # β:    haline contraction (dimensionless, scaled)
    0.035,       # S_0:  reference salinity (psu)
    4.773,       # T_S:  southern temperature
    2.650,       # T_0:  deep ocean temperature
    5.456e6,     # K_N:  North mixing (m³/s)
    5.447e6,     # K_S:  South mixing (m³/s)
    2.79e7,      # λ:    hydraulic constant
    0.39,        # γ:    mixing fraction (southern vs Indo-Pacific)
    5.5e-8,      # μ:    momentum advection
    100*3.15e7,  # Y:    year-to-second scaling × 100
    4.4463e16,   # C:    total salinity content### is this different for 1xCO2?
]

"""
    amoc_params_2xco2

Parameter vector for the 3-box AMOC model at 2xCO2.
"""
const amoc_params_2xco2 = [
    0.0,         # H:    hosing flux (Sv), start at zero
    0.3683e17,   # V_N:  North box volume (m³)
    0.5418e17,   # V_T:  tropical box volume (m³)
    0.6097e17,   # V_S:  southern box volume (m³)
    1.4860e17,   # V_IP: Indo-Pacific box volume (m³)
    9.9250e17,   # V_B:  deep box volume (m³)
    0.034427,    # S_S:  southern salinity (psu, dimensional)
    0.034668,    # S_IP: Indo-Pacific salinity (psu, initial guess; overwritten in ODE)
    0.034538,    # S_B:  deep salinity (psu)
    0.486e6,     # F_N:  freshwater flux North (m³/s)
    -0.997e6,    # F_T:  freshwater flux tropical (m³/s)
    0.12,        # α:    thermal expansion (1/°C)
    790.0,       # β:    haline contraction (dimensionless, scaled)
    0.035,       # S_0:  reference salinity (psu)
    7.919,       # T_S:  southern temperature
    3.870,       # T_0:  deep ocean temperature
    1.762e6,     # K_N:  North mixing (m³/s)
    1.872e6,     # K_S:  South mixing (m³/s)
    1.62e7,      # λ:    hydraulic constant
    0.36,        # γ:    mixing fraction (southern vs Indo-Pacific)
    22e-8,       # μ:    momentum advection
    100*3.15e7,  # Y:    year-to-second scaling × 100
    4.4735e16,   # C:    total salinity content
]

# ─────────────────────────────────────────────────────────────────────────────
# ODE rule
# ─────────────────────────────────────────────────────────────────────────────

"""
    amoc_rule(S, params, t) → SVector{2}

ODE right-hand side for the 3-box AMOC model.

State vector `S = [S_N, S_T]` contains the salinity anomalies (×100, i.e.,
100*(S - S_0)) in the North Atlantic and tropical boxes, respectively.

The overturning flow `q` switches sign at the AMOC collapse. When q ≥ 0 the
circulation is in the "on" (NADW-forming) state; when q < 0 the circulation
is reversed (collapsed state).

Returns dS_N/dt and dS_T/dt.
"""
function amoc_rule(S, params, t)
    S_N, S_T = Tuple(S)
    H, V_N, V_T, V_S, V_IP, V_B, S_S, S_IP, S_B,
        F_N, F_T, α, β, S_0, T_S, T_0, K_N, K_S, λ, γ, μ, Y, C = Tuple(params)

    # Apply hosing: freshwater added to North and tropical boxes
    F_N = F_N + 0.1311e6 * H
    F_T = F_T + 0.6961e6 * H

    # Convert salinities to anomaly representation (×100)
    S_S  = 100 * (S_S - S_0)
    S_B  = 100 * (S_B - S_0)

    # Diagnostic: Indo-Pacific salinity from salt conservation
    S_IP = 100 * (C - (V_N * S_N + V_T * S_T + V_S * S_S + V_B * S_B) / 100 -
                  S_0 * (V_N + V_T + V_S + V_IP + V_B)) / V_IP

    # Overturning strength (m³/s); positive = AMOC on
    q = λ * (α * (T_S - T_0) + β * (S_N - S_S) / 100) / (1 + λ * α * μ)

    if q >= 0
        # AMOC-on state: northward surface flow T→N
        dS_N = Y / V_N * (q * (S_T - S_N) / 100 +
                           K_N * (S_T - S_N) / 100 - F_N * S_0)
        dS_T = Y / V_T * (q * (γ * S_S + (1 - γ) * S_IP - S_T) / 100 +
                           K_S * (S_S - S_T) / 100 +
                           K_N * (S_N - S_T) / 100 - F_T * S_0)
    else
        # AMOC-off state: reversed flow
        dS_N = Y / V_N * (-q * (S_B - S_N) / 100 +
                            K_N * (S_T - S_N) / 100 - F_N * S_0)
        dS_T = Y / V_T * (-q * (S_N - S_T) / 100 +
                            K_S * (S_S - S_T) / 100 +
                            K_N * (S_N - S_T) / 100 - F_T * S_0)
    end

    return SVector(dS_N, dS_T)
end

# ─────────────────────────────────────────────────────────────────────────────
# Diagnostic helper
# ─────────────────────────────────────────────────────────────────────────────

"""
    amoc_strength(S, params) → Float64

Compute the AMOC overturning strength q (m³/s) from a state vector `S = [S_N, S_T]`
and the parameter vector. Positive q corresponds to the AMOC-on state.
"""
function amoc_strength(S, params)
    S_N, S_T = S[1], S[2]
    H, V_N, V_T, V_S, V_IP, V_B, S_S, S_IP, S_B,
        F_N, F_T, α, β, S_0, T_S, T_0, K_N, K_S, λ, γ, μ, Y, C = Tuple(params)

    S_S_anom = 100 * (S_S - S_0)
    q = λ * (α * (T_S - T_0) + β * (S_N - S_S_anom) / 100) / (1 + λ * α * μ)
    return q
end

# ─────────────────────────────────────────────────────────────────────────────
# System constructors
# ─────────────────────────────────────────────────────────────────────────────

"""
    amoc3box_1xco2() → (ds, grid, proximity)

Construct the 3-box AMOC `CoupledODEs` system at 1xCO2.

Returns:
- `ds`: `CoupledODEs` with Rosenbrock23 integrator (stiff solver; the ODE
  has a discontinuity in the Jacobian at q=0 that benefits from stiff methods)
- `grid`: tuple of ranges covering the relevant state space in (S_N, S_T)
- `proximity`: `StrictlyMinimumDistance` instance used with `AttractorsViaProximity`

State space is in salinity anomaly units: S_N, S_T ∈ [-1.5, 1.5] (×100 psu).
Initial condition is near the AMOC-on attractor.
"""
function amoc3box_1xco2()
    params = copy(amoc_params_1xco2)

    # Initial condition: near the AMOC-on state
    u0 = SVector(0.0, 0.4)

    ds = CoupledODEs(amoc_rule, u0, params;
        diffeq = (alg = Rosenbrock23(), abstol = 1e-7, reltol = 1e-7))

    # State-space grid for attractor finding via recurrences.
    # Cell size ≈ 0.01 to match the old ε = 0.01 proximity threshold
    # (in the current Attractors API, grid cell size replaces ε).
    xg = range(-1.0, 1.0; length = 101) #-0.25, 0.15
    yg = range(-0.1,  2.0;  length = 106) #-0.1, 0.7
    grid = (xg, yg)

    # Proximity / mapper keywords (passed as keyword args to AttractorsViaRecurrences)
    # Ttr: worst-case convergence to a fixed point across the grid is ~3400 steps.
    # Use 4000 to ensure all trajectories have settled before recurrence counting.
    proximity = (; ε = 0.01, Ttr = 0.0, Δt = 1.0, stop_at_Δt = true,
                   horizon_limit = 1e2, consecutive_lost_steps = 10_000)

    return ds, grid, proximity
end

"""
    amoc3box_2xco2() → (ds, grid, proximity)

Construct the 3-box AMOC `CoupledODEs` system at 2xCO2 (~2°C warmer).

Returns the same tuple as `amoc3box_1xco2` but with the 2xCO2 parameter set.
Under 2xCO2, the AMOC is generally weaker and the tipping point occurs at a
lower hosing value.
"""
function amoc3box_2xco2()
    params = copy(amoc_params_2xco2)

    u0 = SVector(0.0, 0.4)

    ds = CoupledODEs(amoc_rule, u0, params;
        diffeq = (alg = Rosenbrock23(), abstol = 1e-7, reltol = 1e-7))

    xg = range(-1.0, 1.0; length = 101) #-0.25, 0.15
    yg = range(-0.1,  2.0;  length = 106) #-0.1, 0.7
    grid = (xg, yg)

    # Ttr: worst-case convergence to a fixed point across the grid is ~3400 steps.
    # Use 4000 to ensure all trajectories have settled before recurrence counting.
    proximity = (; ε = 0.01, Ttr = 4000, Δt = 1.0, stop_at_Δt = true,
                   horizon_limit = 1e2, consecutive_lost_steps = 10_000)

    return ds, grid, proximity
end
