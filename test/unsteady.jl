using Test

import CCBlade
import UnsteadyKineticRotorDynamics
import StructArrays

# Verify the filter preserves steady outputs when initialized at steady state.

@testset "unsteady step" begin
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

    Vinf = 10.0
    Omega = 20.0
    rho = 1.0
    ops = [CCBlade.simple_op(Vinf, Omega, r[i], rho) for i in eachindex(r)]

    steady = [CCBlade.solve(rotor, sections[i], ops[i]) for i in eachindex(sections)]

    params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.5, 2.0)
    state = UnsteadyKineticRotorDynamics.UnsteadyState(sections, ops; V_wake_old = Vinf)
    state.outputs = StructArrays.StructArray(steady)

    out = UnsteadyKineticRotorDynamics.unsteady_step!(state, rotor, sections, ops, params; dt = 0.1)
    loads = UnsteadyKineticRotorDynamics.rotor_loads(rotor, sections, out; omega = Omega)
    T, Q = CCBlade.thrusttorque(rotor, sections, out)

    @test isapprox(state.time, 0.1; atol = 1e-12)
    steady_Np = [steady[i].Np for i in eachindex(steady)]
    steady_Tp = [steady[i].Tp for i in eachindex(steady)]

    @test all(isapprox.(out.Np, steady_Np; atol = 1e-12))
    @test all(isapprox.(out.Tp, steady_Tp; atol = 1e-12))
    @test isapprox(loads.thrust_n, T; atol = 1e-12)
    @test isapprox(loads.torque_nm, Q; atol = 1e-12)
    @test isapprox(loads.shaft_power_w, Q * Omega; atol = 1e-12)
end

@testset "simple rotor construction" begin
    rotor_components = UnsteadyKineticRotorDynamics.simple_blade_element_rotor(
        rotor_radius = 1.0,
        hub_radius = 0.1,
        blades = 2,
        n_sections = 3,
    )

    @test length(rotor_components.radii) == 3
    @test length(rotor_components.sections) == 3
    @test rotor_components.rotor.Rtip == 1.0
end
