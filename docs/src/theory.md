# Theory

UnsteadyKineticRotorDynamics computes steady CCBlade outputs at each time step,
then applies a low-pass filter to approximate unsteady kinetic rotor response.
Two time constants are used:

- `tau_near`: filters sectional outputs (loads, induction factors, velocities).
- `tau_far`: filters the wake speed used to scale the time constants at the next step.

## Low-Pass Filter

For any sectional quantity `q`, the filtered update is

```
q_k = w * q_{k-1} + (1 - w) * q_steady
w = exp(-dt / tau_near)
```

The time constants are nondimensional and scaled by the rotor tip radius and the previous wake speed:

```
tau_near_dim = tau_near * Rtip / V_wake_old
tau_far_dim  = tau_far  * Rtip / V_wake_old
```

## Wake Speed Update

The wake speed is updated using a similar low-pass form:

```
V_wake_new = V_wake_old * w_far + Vx_mean * (1 - 2 * a_mean) * (1 - w_far)
```

where `a_mean` is the mean axial induction factor (sign-adjusted for turbine/propeller convention), and `w_far = exp(-dt / tau_far_dim)`.

## Base-Motion Inflow

The optional pole-motion helper computes the hub velocity from base translation and rotation about an arm. The relative wind at the hub is:

```
Vrel = Vhub - Vhub_vel
```

Yaw and tilt are then updated using the relative wind direction, and roll is applied as an azimuth shift. This prevents double counting of velocity changes while still accounting for base motion in the inflow geometry.
