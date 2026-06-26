using Test

import CCBlade
import UnsteadyKineticRotorDynamics

# Sine-modulated inflow to exercise unsteady stepping.

@testset "unsteady sine" begin
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

    rho = 1.0
    Omega = 20.0
    Vinf_nominal = 10.0

    params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.5, 2.0)
    ops0 = [CCBlade.simple_op(Vinf_nominal, Omega, r[i], rho) for i in eachindex(r)]
    state = UnsteadyKineticRotorDynamics.UnsteadyState(sections, ops0; V_wake_old = Vinf_nominal)

    nsteps = 100
    dt = 0.1

    thrust = zeros(nsteps)

    for i in 1:nsteps
        Vinf = Vinf_nominal * sin(pi * (i - 1) / (nsteps - 1))
        ops = [CCBlade.simple_op(Vinf, Omega, r[i], rho) for i in eachindex(r)]

        snapshot = UnsteadyKineticRotorDynamics.unsteady_loads_step!(
            state, rotor, sections, ops, params; dt = dt, omega = Omega)
        thrust[i] = snapshot.thrust_n
    end

    @test isapprox(state.time, nsteps * dt; atol = 1e-12)
    @test all(isfinite.(thrust))
    @test any(abs.(thrust) .> 0)
end
