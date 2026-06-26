__precompile__()

module UnsteadyKineticRotorDynamics

import CCBlade
import StructArrays
using StructArrays: StructArray

export UnsteadyParams, UnsteadyState
export simple_blade_element_rotor, windturbine_op_motion
export rotor_loads, unsteady_loads_step, unsteady_loads_step!, unsteady_step, unsteady_step!

# Low-pass unsteady stepping utilities layered on steady CCBlade solves.

const _ifw_loaded = Ref(false)

@inline function _mean(x)
    return sum(x) / length(x)
end

@inline function _op_field(ops, field::Symbol)
    if ops isa StructArrays.StructArray
        return StructArrays.component(ops, field)
    end
    return [getfield(op, field) for op in ops]
end

@inline function _section_field(sections, field::Symbol)
    if sections isa StructArrays.StructArray
        return StructArrays.component(sections, field)
    end
    return [getfield(section, field) for section in sections]
end

@inline function _output_field(outputs, field::Symbol)
    if outputs isa StructArrays.StructArray
        return StructArrays.component(outputs, field)
    end
    return [getfield(output, field) for output in outputs]
end

@inline function _as_tuple3(v)
    return (v[1], v[2], v[3])
end

@inline function _cross(a, b)
    return (a[2] * b[3] - a[3] * b[2],
        a[3] * b[1] - a[1] * b[3],
        a[1] * b[2] - a[2] * b[1])
end

@inline function _rotate_rpy(v, angles)
    x, y, z = _as_tuple3(v)
    roll, pitch, yaw = _as_tuple3(angles)

    cr = cos(roll)
    sr = sin(roll)
    cp = cos(pitch)
    sp = sin(pitch)
    cy = cos(yaw)
    sy = sin(yaw)

    xw = cy * cp * x + (-sy * cr + cy * sp * sr) * y + (sy * sr + cy * sp * cr) * z
    yw = sy * cp * x + (cy * cr + sy * sp * sr) * y + (-cy * sr + sy * sp * cr) * z
    zw = -sp * x + cp * sr * y + cp * cr * z

    return (xw, yw, zw)
end

@inline function _zero_outputs(zero_val)
    return CCBlade.Outputs(zero_val, zero_val, zero_val, zero_val, zero_val,
        zero_val, zero_val, zero_val, zero_val, zero_val,
        zero_val, zero_val, zero_val, zero_val, zero_val)
end

function _init_outputs(n::Integer, init_output::CCBlade.Outputs)
    return StructArray([init_output for _ in 1:n])
end

function _ifw_module()
    if !_ifw_loaded[]
        try
            @eval import OWENSOpenFASTWrappers
        catch err
            error("OWENSOpenFASTWrappers is required when ifw=true. Add it to your environment and try again.")
        end
        _ifw_loaded[] = true
    end
    return OWENSOpenFASTWrappers
end

function _default_ifw_positions(sections, rotor, azimuth, ifw_center, ifw_z, rotation)
    r = _section_field(sections, :r)
    n = length(r)
    R = rotor.Rtip
    if iszero(R)
        error("rotor.Rtip must be nonzero for ifw positioning.")
    end

    if isa(azimuth, Number)
        sa = sin(azimuth)
        ca = cos(azimuth)
        x = rotation .* sa .* (r ./ R)
        y = rotation .* ca .* (r ./ R)
    else
        if length(azimuth) != n
            error("azimuth must be a scalar or match the number of sections.")
        end
        sa = sin.(azimuth)
        ca = cos.(azimuth)
        x = rotation .* sa .* (r ./ R)
        y = rotation .* ca .* (r ./ R)
    end

    return [(ifw_center[1] + x[i], ifw_center[2] + y[i], ifw_z) for i in 1:n]
end

function _ifw_ops(ops, sections, rotor, time;
        azimuth = 0.0,
        ifw_positions = nothing,
        ifw_center = (0.0, 0.0),
        ifw_z = 0.0,
        ifw_rotation = nothing,
        ifw_mode = :axial)
    n = length(sections)
    if length(ops) != n
        error("sections and ops must have the same length.")
    end

    ifwmod = _ifw_module()

    if isnothing(ifw_rotation)
        mean_vy = _mean(_op_field(ops, :Vy))
        rotation = sign(mean_vy)
        if iszero(rotation)
            rotation = one(mean_vy)
        end
    else
        rotation = ifw_rotation
    end
    positions = ifw_positions === nothing ?
                _default_ifw_positions(
        sections, rotor, azimuth, ifw_center, ifw_z, rotation) :
                ifw_positions

    if length(positions) != n
        error("ifw_positions must match the number of sections.")
    end

    Vx_base = _op_field(ops, :Vx)
    Vy_base = _op_field(ops, :Vy)
    Vx = similar(Vx_base)
    Vy = similar(Vy_base)

    for i in 1:n
        pos = positions[i]
        vel = ifwmod.ifwcalcoutput([pos[1], pos[2], pos[3]], time)
        vx = vel[1]
        vy = vel[2]

        if ifw_mode == :replace
            Vx[i] = vx
            Vy[i] = vy
        elseif ifw_mode == :add
            Vx[i] = Vx_base[i] + vx
            Vy[i] = Vy_base[i] + vy
        elseif ifw_mode == :axial
            Vx[i] = vx
            Vy[i] = Vy_base[i]
        else
            error("ifw_mode must be :axial, :add, or :replace.")
        end
    end

    rho = _op_field(ops, :rho)
    pitch = _op_field(ops, :pitch)
    mu = _op_field(ops, :mu)
    asound = _op_field(ops, :asound)

    return [CCBlade.OperatingPoint(Vx[i], Vy[i], rho[i], pitch[i], mu[i], asound[i])
            for i in 1:n]
end

function _filter_outputs!(prev, cur, w)
    wp = one(w) - w
    prev_Np = _output_field(prev, :Np)
    prev_Tp = _output_field(prev, :Tp)
    prev_a = _output_field(prev, :a)
    prev_ap = _output_field(prev, :ap)
    prev_u = _output_field(prev, :u)
    prev_v = _output_field(prev, :v)
    prev_phi = _output_field(prev, :phi)
    prev_alpha = _output_field(prev, :alpha)
    prev_W = _output_field(prev, :W)
    prev_cl = _output_field(prev, :cl)
    prev_cd = _output_field(prev, :cd)
    prev_cn = _output_field(prev, :cn)
    prev_ct = _output_field(prev, :ct)
    prev_F = _output_field(prev, :F)
    prev_G = _output_field(prev, :G)

    cur_Np = _output_field(cur, :Np)
    cur_Tp = _output_field(cur, :Tp)
    cur_a = _output_field(cur, :a)
    cur_ap = _output_field(cur, :ap)
    cur_u = _output_field(cur, :u)
    cur_v = _output_field(cur, :v)
    cur_phi = _output_field(cur, :phi)
    cur_alpha = _output_field(cur, :alpha)
    cur_W = _output_field(cur, :W)
    cur_cl = _output_field(cur, :cl)
    cur_cd = _output_field(cur, :cd)
    cur_cn = _output_field(cur, :cn)
    cur_ct = _output_field(cur, :ct)
    cur_F = _output_field(cur, :F)
    cur_G = _output_field(cur, :G)

    prev_Np .= prev_Np .* w .+ cur_Np .* wp
    prev_Tp .= prev_Tp .* w .+ cur_Tp .* wp
    prev_a .= prev_a .* w .+ cur_a .* wp
    prev_ap .= prev_ap .* w .+ cur_ap .* wp
    prev_u .= prev_u .* w .+ cur_u .* wp
    prev_v .= prev_v .* w .+ cur_v .* wp
    prev_phi .= prev_phi .* w .+ cur_phi .* wp
    prev_alpha .= prev_alpha .* w .+ cur_alpha .* wp
    prev_W .= prev_W .* w .+ cur_W .* wp
    prev_cl .= prev_cl .* w .+ cur_cl .* wp
    prev_cd .= prev_cd .* w .+ cur_cd .* wp
    prev_cn .= prev_cn .* w .+ cur_cn .* wp
    prev_ct .= prev_ct .* w .+ cur_ct .* wp
    prev_F .= prev_F .* w .+ cur_F .* wp
    prev_G .= prev_G .* w .+ cur_G .* wp
    return prev
end

"""
    UnsteadyParams(tau_near, tau_far; ifw=false)
    UnsteadyParams(tau; ifw=false)

Parameters for the unsteady low-pass filter.
`tau_near` and `tau_far` are nondimensional time constants scaled by `Rtip / V_wake_old` each step.

# Examples
```julia
import UnsteadyKineticRotorDynamics

params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.3, 3.0)
params_ifw = UnsteadyKineticRotorDynamics.UnsteadyParams((0.3, 3.0); ifw=true)
```
"""
struct UnsteadyParams{TF, TB}
    tau_near::TF
    tau_far::TF
    ifw::TB
end

function UnsteadyParams(tau_near, tau_far; ifw = false)
    if !(tau_near > 0 && tau_far > 0)
        error("tau_near and tau_far must be positive.")
    end
    return UnsteadyParams(promote(tau_near, tau_far)..., ifw)
end

function UnsteadyParams(tau::AbstractVector; ifw = false)
    if length(tau) != 2
        error("tau must have length 2.")
    end
    return UnsteadyParams(tau[1], tau[2]; ifw = ifw)
end

UnsteadyParams(tau::NTuple{2}; ifw = false) = UnsteadyParams(tau[1], tau[2]; ifw = ifw)

"""
    simple_blade_element_rotor(; kwargs...)

Create a compact blade-element rotor and section set using the package's CCBlade
backend. This keeps downstream packages from constructing backend-specific rotor
objects directly.
"""
function simple_blade_element_rotor(;
        rotor_radius,
        hub_radius = 0.1 * rotor_radius,
        blades::Integer = 3,
        n_sections::Integer = 5,
        chord_fraction = 0.08,
        theta_rad = 0.0,
        lift_slope = 6.2,
        drag_coefficient = 0.01,
        precone_rad = 0.0,
        turbine::Bool = true)
    if n_sections < 2
        throw(ArgumentError("n_sections must be at least 2"))
    end
    if !(rotor_radius > hub_radius > zero(rotor_radius))
        throw(ArgumentError("rotor_radius must be greater than positive hub_radius"))
    end

    radii = collect(range(hub_radius + (rotor_radius - hub_radius) / (2n_sections),
        rotor_radius - (rotor_radius - hub_radius) / (2n_sections), length = n_sections))
    chord = chord_fraction * rotor_radius

    function affun(alpha, Re, M)
        return lift_slope * alpha, drag_coefficient
    end

    rotor = CCBlade.Rotor(hub_radius, rotor_radius, blades; precone = precone_rad, turbine = turbine)
    sections = [CCBlade.Section(r, chord, theta_rad, affun) for r in radii]
    return (rotor = rotor, sections = sections, radii = radii)
end

"""
    UnsteadyState(sections, ops; init_output=nothing, V_wake_old=nothing, time=nothing)
    UnsteadyState(n; init_output=Outputs(), V_wake_old=1.0, time=0.0)

State for unsteady stepping. `outputs` stores filtered CCBlade outputs.

# Examples
```julia
import CCBlade
import UnsteadyKineticRotorDynamics

r = [0.2, 0.6, 0.9]
chord = fill(0.1, 3)
theta = fill(0.0, 3)

affun = (alpha, Re, M) -> (6.2 * alpha, 0.01)
sections = [CCBlade.Section(r[i], chord[i], theta[i], affun) for i in eachindex(r)]
ops = [CCBlade.simple_op(10.0, 20.0, r[i], 1.0) for i in eachindex(r)]

state = UnsteadyKineticRotorDynamics.UnsteadyState(sections, ops; V_wake_old=10.0)
state2 = UnsteadyKineticRotorDynamics.UnsteadyState(length(sections))
```
"""
mutable struct UnsteadyState{TO, TF}
    outputs::TO
    V_wake_old::TF
    time::TF
end

function UnsteadyState(
        n::Integer; init_output = CCBlade.Outputs(), V_wake_old = 1.0, time = 0.0)
    outputs = _init_outputs(n, init_output)
    return UnsteadyState(outputs, V_wake_old, time)
end

function UnsteadyState(sections::AbstractVector, ops::AbstractVector;
        init_output = nothing,
        V_wake_old = nothing,
        time = nothing)
    n = length(sections)
    Vx = _op_field(ops, :Vx)
    if init_output === nothing
        init_output = _zero_outputs(zero(eltype(Vx)))
    end
    outputs = _init_outputs(n, init_output)
    V_wake_val = isnothing(V_wake_old) ? _mean(Vx) : V_wake_old
    time_val = isnothing(time) ? zero(V_wake_val) : time
    return UnsteadyState(outputs, V_wake_val, time_val)
end

"""
    windturbine_op_motion(Vhub, Omega, pitch, r, precone, yaw, tilt, azimuth,
        hubHt, shearExp, rho; base_pos=(0,0,0), base_vel=(0,0,0),
        base_angles=(0,0,0), base_omega=(0,0,0), arm=(0,0,0),
        mu=one(rho), asound=one(rho))

Compute windturbine operating points with hub motion from a 6-DOF base motion
applied about the base of an arm (pole). The arm vector is defined in the base
frame and represents the nominal hub location. `base_pos` and `base_vel` are
displacements and velocities of the base in the wind frame (x along wind,
y lateral, z up). `base_angles` are roll, pitch, yaw (rad) about x, y, z using
a Z-Y-X rotation, and `base_omega` is the angular velocity vector in the wind
frame (rad/s). `base_omega` is not generally equal to the time derivatives of
`base_angles` except in the small-angle limit.

The hub velocity is subtracted from the inflow to avoid double counting.
Yaw and tilt are updated by the base angles and by the relative wind direction.
Returns `(ops, info)` where `info` includes `hub_pos`, `hub_vel`, `yaw_eff`, `tilt_eff`,
`azimuth_eff`, `Vhub_eff`, and `hubHt_eff`.

# Examples
```julia
import CCBlade
import UnsteadyKineticRotorDynamics

r = [5.0, 10.0, 15.0]

ops, info = UnsteadyKineticRotorDynamics.windturbine_op_motion(
    10.0, 2.0, 0.0, r, 0.0, 0.0, 0.0, 0.0, 100.0, 0.2, 1.225;
    base_pos=(0.0, 0.0, 1.0),
    base_vel=(0.0, 0.0, -0.1),
    base_angles=(0.0, 0.0, 0.0),
    base_omega=(0.0, 0.0, 0.0),
    arm=(0.0, 0.0, 100.0),
)
```
"""
function windturbine_op_motion(Vhub, Omega, pitch, r, precone, yaw, tilt, azimuth,
        hubHt, shearExp, rho; base_pos = (0.0, 0.0, 0.0),
        base_vel = (0.0, 0.0, 0.0),
        base_angles = (0.0, 0.0, 0.0),
        base_omega = (0.0, 0.0, 0.0),
        arm = (0.0, 0.0, 0.0),
        mu = one(rho),
        asound = one(rho))
    base_pos_t = _as_tuple3(base_pos)
    base_vel_t = _as_tuple3(base_vel)
    base_angles_t = _as_tuple3(base_angles)
    base_omega_t = _as_tuple3(base_omega)
    arm_t = _as_tuple3(arm)

    arm_rot = _rotate_rpy(arm_t, base_angles_t)
    hub_pos = (base_pos_t[1] + arm_rot[1],
        base_pos_t[2] + arm_rot[2],
        base_pos_t[3] + arm_rot[3])

    hub_vel = (base_vel_t[1], base_vel_t[2], base_vel_t[3])
    omega_cross_arm = _cross(base_omega_t, arm_rot)
    hub_vel = (hub_vel[1] + omega_cross_arm[1],
        hub_vel[2] + omega_cross_arm[2],
        hub_vel[3] + omega_cross_arm[3])

    hub_disp_z = base_pos_t[3] + arm_rot[3] - arm_t[3]
    hubHt_eff = hubHt + hub_disp_z

    # Relative wind in the hub frame (wind minus hub velocity).
    Vrel_x = Vhub - hub_vel[1]
    Vrel_y = -hub_vel[2]
    Vrel_z = -hub_vel[3]
    Vhub_eff = sqrt(Vrel_x * Vrel_x + Vrel_y * Vrel_y + Vrel_z * Vrel_z)

    wind_yaw = atan(Vrel_y, Vrel_x)
    wind_tilt = atan(Vrel_z, hypot(Vrel_x, Vrel_y))

    roll = base_angles_t[1]
    pitch_b = base_angles_t[2]
    yaw_b = base_angles_t[3]

    yaw_eff = yaw + yaw_b - wind_yaw
    tilt_eff = tilt + pitch_b - wind_tilt
    # Roll about the wind axis shifts the azimuth reference.
    azimuth_eff = azimuth + roll

    if isa(r, AbstractVector)
        ops = [CCBlade.windturbine_op(Vhub_eff, Omega, pitch, r[i], precone,
            yaw_eff, tilt_eff, azimuth_eff, hubHt_eff, shearExp, rho, mu, asound)
            for i in eachindex(r)]
    else
        ops = CCBlade.windturbine_op(Vhub_eff, Omega, pitch, r, precone, yaw_eff,
            tilt_eff, azimuth_eff, hubHt_eff, shearExp, rho, mu, asound)
    end

    info = (hub_pos = hub_pos,
        hub_vel = hub_vel,
        yaw_eff = yaw_eff,
        tilt_eff = tilt_eff,
        azimuth_eff = azimuth_eff,
        Vhub_eff = Vhub_eff,
        hubHt_eff = hubHt_eff)

    return ops, info
end

"""
    unsteady_step!(state, rotor, sections, ops, params;
        dt,
        azimuth=0.0,
        ifw=params.ifw,
        ifw_positions=nothing,
        ifw_center=(0.0, 0.0),
        ifw_z=0.0,
        ifw_rotation=nothing,
        ifw_mode=:axial,
        time=nothing)

Advance the unsteady model by one time step using a low-pass filter on the steady CCBlade outputs.
If `ifw` is true, inflow velocities are updated using `OWENSOpenFASTWrappers.ifwcalcoutput`.
Default `ifw_positions` uses rotor-plane coordinates scaled by `r / Rtip` (override with `ifw_positions`).
`ifw_mode` controls how inflow wind is applied:
- `:axial` replaces only `Vx` (keeps `Vy` from `ops`).
- `:add` adds inflow components to `ops` (treats them as perturbations).
- `:replace` replaces both `Vx` and `Vy`.

Returns the filtered outputs and mutates `state`.

# Examples
```julia
import CCBlade
import UnsteadyKineticRotorDynamics

Rhub = 0.1
Rtip = 1.0
rotor = CCBlade.Rotor(Rhub, Rtip, 2)

r = [0.2, 0.6, 0.9]
chord = fill(0.1, 3)
theta = fill(0.0, 3)

affun = (alpha, Re, M) -> (6.2 * alpha, 0.01)
sections = [CCBlade.Section(r[i], chord[i], theta[i], affun) for i in eachindex(r)]
ops = [CCBlade.simple_op(10.0, 20.0, r[i], 1.0) for i in eachindex(r)]

params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.3, 3.0)
state = UnsteadyKineticRotorDynamics.UnsteadyState(sections, ops; V_wake_old=10.0)
out = UnsteadyKineticRotorDynamics.unsteady_step!(state, rotor, sections, ops, params; dt=0.1)
```
"""
function unsteady_step!(state::UnsteadyState, rotor, sections, ops, params::UnsteadyParams;
        dt,
        azimuth = 0.0,
        ifw = params.ifw,
        ifw_positions = nothing,
        ifw_center = (0.0, 0.0),
        ifw_z = 0.0,
        ifw_rotation = nothing,
        ifw_mode = :axial,
        time = nothing)
    n = length(sections)
    if length(ops) != n
        error("sections and ops must have the same length.")
    end
    if length(state.outputs) != n
        error("state.outputs length must match the number of sections.")
    end

    if !(dt > 0)
        error("dt must be positive.")
    end
    if !(rotor.Rtip > 0)
        error("rotor.Rtip must be positive.")
    end

    t = isnothing(time) ? state.time + dt : time

    ops_use = ifw ?
              _ifw_ops(ops, sections, rotor, t;
        azimuth = azimuth,
        ifw_positions = ifw_positions,
        ifw_center = ifw_center,
        ifw_z = ifw_z,
        ifw_rotation = ifw_rotation,
        ifw_mode = ifw_mode) :
              ops

    steady_out = [CCBlade.solve(rotor, sections[i], ops_use[i])
                  for i in eachindex(sections)]

    a_arr = _output_field(steady_out, :a)
    a_phys = rotor.turbine ? -a_arr : a_arr
    a_mean = _mean(a_phys)
    Vx_mean = _mean(_op_field(ops_use, :Vx))

    V_wake_prev = state.V_wake_old
    V_wake_mag = abs(V_wake_prev)
    small = oftype(V_wake_mag, 1e-6)
    V_wake_ref = ifelse(V_wake_mag > small, V_wake_mag, small)

    tau_near = params.tau_near * rotor.Rtip / V_wake_ref
    tau_far = params.tau_far * rotor.Rtip / V_wake_ref

    w_near = exp(-dt / tau_near)
    w_far = exp(-dt / tau_far)

    _filter_outputs!(state.outputs, steady_out, w_near)

    state.V_wake_old = V_wake_prev * w_far +
                       Vx_mean * (1 - 2 * a_mean) * (one(w_far) - w_far)
    state.time = t

    return state.outputs
end

"""
    unsteady_step(rotor, sections, ops, params; dt, kwargs...)

Convenience wrapper that creates a new `UnsteadyState` and returns `(outputs, state)`.

# Examples
```julia
import CCBlade
import UnsteadyKineticRotorDynamics

Rhub = 0.1
Rtip = 1.0
rotor = CCBlade.Rotor(Rhub, Rtip, 2)

r = [0.2, 0.6, 0.9]
chord = fill(0.1, 3)
theta = fill(0.0, 3)

affun = (alpha, Re, M) -> (6.2 * alpha, 0.01)
sections = [CCBlade.Section(r[i], chord[i], theta[i], affun) for i in eachindex(r)]
ops = [CCBlade.simple_op(10.0, 20.0, r[i], 1.0) for i in eachindex(r)]

params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.3, 3.0)
out, state = UnsteadyKineticRotorDynamics.unsteady_step(rotor, sections, ops, params; dt=0.1)
```
"""
function unsteady_step(rotor, sections, ops, params::UnsteadyParams; dt, kwargs...)
    state = UnsteadyState(sections, ops)
    outputs = unsteady_step!(state, rotor, sections, ops, params; dt = dt, kwargs...)
    return outputs, state
end

"""
    rotor_loads(rotor, sections, outputs; omega)

Integrate section outputs into rotor thrust, torque, and shaft power. The sign
convention matches the CCBlade backend; downstream systems can cap or rectify
shaft power according to their generator and control model.
"""
function rotor_loads(rotor, sections, outputs; omega)
    thrust_n, torque_nm = CCBlade.thrusttorque(rotor, sections, outputs)
    shaft_power_w = torque_nm * omega
    return (thrust_n = thrust_n, torque_nm = torque_nm, shaft_power_w = shaft_power_w)
end

"""
    unsteady_loads_step!(state, rotor, sections, ops, params; dt, omega, kwargs...)

Advance the unsteady rotor state and return integrated rotor loads and filtered
outputs. This is the preferred boundary for system-level models that need
thrust, torque, and shaft power without depending on the CCBlade backend.
"""
function unsteady_loads_step!(state::UnsteadyState, rotor, sections, ops, params::UnsteadyParams;
        dt, omega, kwargs...)
    outputs = unsteady_step!(state, rotor, sections, ops, params; dt = dt, kwargs...)
    loads = rotor_loads(rotor, sections, outputs; omega = omega)
    return (
        outputs = outputs,
        thrust_n = loads.thrust_n,
        torque_nm = loads.torque_nm,
        shaft_power_w = loads.shaft_power_w,
    )
end

"""
    unsteady_loads_step(rotor, sections, ops, params; dt, omega, kwargs...)

Convenience wrapper that creates a new `UnsteadyState` and returns
`(snapshot, state)`, where `snapshot` includes integrated rotor loads.
"""
function unsteady_loads_step(rotor, sections, ops, params::UnsteadyParams; dt, omega, kwargs...)
    state = UnsteadyState(sections, ops)
    snapshot = unsteady_loads_step!(state, rotor, sections, ops, params;
        dt = dt, omega = omega, kwargs...)
    return snapshot, state
end

end
