using Test

import UnsteadyKineticRotorDynamics

@testset "wind turbine helper utilities" begin
    rotor_parts = UnsteadyKineticRotorDynamics.simple_blade_element_rotor(
        rotor_radius = 1.0,
        hub_radius = 0.1,
        blades = 3,
        n_sections = 4,
        chord_fraction = 0.1,
    )
    params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.3, 3.0)

    Vhub = 6.0
    tsr = 4.0
    pitch = 0.0
    precone = 0.0
    yaw = 0.0
    tilt = 0.0
    hub_height = 10.0
    shear_exp = 0.0
    rho = 1.225
    dt = 0.05

    ops, info, omega = UnsteadyKineticRotorDynamics.windturbine_op_motion_tsr(
        Vhub,
        tsr,
        pitch,
        rotor_parts.radii,
        precone,
        yaw,
        tilt,
        0.0,
        hub_height,
        shear_exp,
        rho;
        rotor_radius = rotor_parts.rotor.Rtip,
        arm = (0.0, 0.0, hub_height),
    )

    @test length(ops) == length(rotor_parts.sections)
    @test isapprox(omega, tsr * info.Vhub_eff / rotor_parts.rotor.Rtip; atol = 1e-12)

    loads = UnsteadyKineticRotorDynamics.settled_windturbine_loads(
        rotor_parts.rotor,
        rotor_parts.sections,
        rotor_parts.radii;
        tsr = tsr,
        Vhub = Vhub,
        hub_height = hub_height,
        dt = dt,
        params = params,
        settling_steps = 3,
        rho = rho,
    )
    @test isfinite(loads.shaft_power_w)
    @test isfinite(loads.thrust_n)

    calibration = UnsteadyKineticRotorDynamics.calibrate_windturbine_tsr(
        rotor_parts.rotor,
        rotor_parts.sections,
        rotor_parts.radii;
        tsr_values = 3.0:0.5:5.0,
        Vhub = Vhub,
        hub_height = hub_height,
        dt = dt,
        params = params,
        settling_steps = 3,
        rho = rho,
        target_power_w = 100.0,
    )
    @test calibration.optimal_tsr in collect(3.0:0.5:5.0)
    @test isfinite(calibration.steady_raw_power_w)
    @test 0.0 < calibration.power_scale <= 1.0

    azimuths = (0.0, 2pi / 3, 4pi / 3)
    operating_point = azimuth -> UnsteadyKineticRotorDynamics.windturbine_op_motion_tsr(
        Vhub,
        tsr,
        pitch,
        rotor_parts.radii,
        precone,
        yaw,
        tilt,
        azimuth,
        hub_height,
        shear_exp,
        rho;
        rotor_radius = rotor_parts.rotor.Rtip,
        arm = (0.0, 0.0, hub_height),
    )
    states = UnsteadyKineticRotorDynamics.azimuth_unsteady_states(
        rotor_parts.sections, azimuths, operating_point)
    avg_loads = UnsteadyKineticRotorDynamics.azimuth_averaged_unsteady_loads_step!(
        states,
        rotor_parts.rotor,
        rotor_parts.sections,
        params,
        azimuths,
        operating_point;
        dt = dt,
        shaft_power_scale = 0.5,
    )

    @test length(states) == length(azimuths)
    @test isfinite(avg_loads.raw_shaft_power_w)
    @test isfinite(avg_loads.thrust_n)
    @test isapprox(avg_loads.shaft_power_w, 0.5 * avg_loads.raw_shaft_power_w;
        atol = 1e-12)
end
