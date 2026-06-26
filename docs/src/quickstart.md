# Quickstart

UnsteadyKineticRotorDynamics wraps steady CCBlade solves with a low-pass time filter. This quickstart shows a minimal setup with a simple rotor and three blade sections.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/kevmoor/UnsteadyKineticRotorDynamics.jl")
```

## Basic Unsteady Stepping

```julia
import CCBlade
import UnsteadyKineticRotorDynamics

# Rotor and sections
Rhub = 0.1
Rtip = 1.0
B = 2
rotor = CCBlade.Rotor(Rhub, Rtip, B)

r = [0.2, 0.6, 0.9]
chord = fill(0.1, 3)
theta = fill(0.0, 3)

function affun(alpha, Re, M)
    return 6.2 * alpha, 0.01
end

sections = [CCBlade.Section(r[i], chord[i], theta[i], affun) for i in eachindex(r)]

# Operating conditions
rho = 1.0
Omega = 20.0
Vinf = 10.0
ops = [CCBlade.simple_op(Vinf, Omega, r[i], rho) for i in eachindex(r)]

# Unsteady setup
params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.3, 3.0)
state = UnsteadyKineticRotorDynamics.UnsteadyState(sections, ops; V_wake_old=Vinf)

dt = 0.1
for k in 1:50
    out = UnsteadyKineticRotorDynamics.unsteady_step!(state, rotor, sections, ops, params; dt=dt)
end
```

## Base Motion (Pole) Inflow

```julia
ops, info = UnsteadyKineticRotorDynamics.windturbine_op_motion(
    10.0, 20.0, 0.0, r, 0.0, 0.0, 0.0, 0.0, 10.0, 0.2, 1.225;
    base_pos=(0.0, 0.0, 0.5),
    base_vel=(0.0, 0.0, -0.1),
    base_angles=(0.0, 0.0, 0.0),
    base_omega=(0.0, 0.0, 0.0),
    arm=(0.0, 0.0, 10.0),
)
```

The returned `info` named tuple includes `hub_pos`, `hub_vel`, `yaw_eff`, `tilt_eff`, `azimuth_eff`, `Vhub_eff`, and `hubHt_eff`.

## Hydrokinetic Use

The same rotor wrapper can be used for marine hydrokinetic studies when the fluid
properties and geometry convention are set consistently. Use water density for
`rho`, water-current speed for inflow, and orient the tower/arm convention so the
hub location is below the floating platform rather than above it. SIRENOpt uses
this boundary to share rotor inflow, speed, torque, power, and load quantities
between wind and hydrokinetic placeholder models.
