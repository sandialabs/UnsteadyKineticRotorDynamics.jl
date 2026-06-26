using Test

import CCBlade
import UnsteadyKineticRotorDynamics
import ForwardDiff
import FiniteDiff

# Set true to see plots when running interactively.
makeplots = true

@testset "pole motion" begin
    # --- rotor definition ---
    Rhub = 0.1
    Rtip = 1.0
    B = 2
    precone = 0.0
    rotor = CCBlade.Rotor(Rhub, Rtip, B, precone=precone, turbine=true)

    # --- section definitions ---
    r = [0.2, 0.6, 0.9]
    chord = fill(0.1, 3)
    theta = fill(0.0, 3)

    function affun(alpha, Re, M)
        return 6.2 * alpha, 0.01
    end

    sections = [CCBlade.Section(r[i], chord[i], theta[i], affun) for i in eachindex(r)]

    # --- operating conditions ---
    Vhub = 10.0
    Omega = 20.0
    pitch = 0.0
    yaw = 0.0
    tilt = 0.0
    azimuth = 0.0
    hubHt = 10.0
    shearExp = 0.0
    rho = 1.0

    arm = (0.0, 0.0, hubHt)

    params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.3, 3.0)
    ops0, _ = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        yaw, tilt, azimuth, hubHt, shearExp, rho; arm=arm)
    state = UnsteadyKineticRotorDynamics.UnsteadyState(sections, ops0; V_wake_old=Vhub)

    nsteps = 50
    dt = 0.1
    omega_m = 2.0 * pi / (nsteps * dt)

    thrust = zeros(nsteps)
    vhub_eff = zeros(nsteps)
    power = zeros(nsteps)

    for i in 1:nsteps
        t = (i - 1) * dt
        heave = 0.2 * sin(omega_m * t)
        heave_dot = 0.2 * omega_m * cos(omega_m * t)

        roll = 1.0 * pi / 180 * sin(omega_m * t)
        pitch_b = 2.0 * pi / 180 * sin(omega_m * t)
        yaw_b = 0.5 * pi / 180 * sin(omega_m * t)

        roll_dot = 1.0 * pi / 180 * omega_m * cos(omega_m * t)
        pitch_dot = 2.0 * pi / 180 * omega_m * cos(omega_m * t)
        yaw_dot = 0.5 * pi / 180 * omega_m * cos(omega_m * t)

        base_pos = (0.0, 0.0, heave)
        base_vel = (0.0, 0.0, heave_dot)
        base_angles = (roll, pitch_b, yaw_b)
        base_omega = (roll_dot, pitch_dot, yaw_dot)

        ops, info = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
            yaw, tilt, azimuth, hubHt, shearExp, rho;
            base_pos=base_pos,
            base_vel=base_vel,
            base_angles=base_angles,
            base_omega=base_omega,
            arm=arm)

        out = UnsteadyKineticRotorDynamics.unsteady_step!(state, rotor, sections, ops, params; dt=dt)
        T, Q = CCBlade.thrusttorque(rotor, sections, out)
        thrust[i] = T
        vhub_eff[i] = info.Vhub_eff
        power[i] = Q * Omega
    end

    @test all(isfinite.(thrust))
    @test all(isfinite.(vhub_eff))
    @test all(isfinite.(power))
    @test any(abs.(thrust) .> 0)

    if makeplots == true
        try
            @eval import Plots
            figdir = joinpath(@__DIR__, "..", "figs")
            mkpath(figdir)
            t = (0:nsteps-1) .* dt
            p1 = Plots.plot(t, vhub_eff, xlabel="Time (s)", ylabel="Vhub eff (m/s)",
                label="Vhub eff", background_color=:transparent,
                background_color_inside=:transparent, background_color_outside=:transparent)
            p2 = Plots.plot(t, thrust, xlabel="Time (s)", ylabel="Thrust (N)",
                label="Thrust", background_color=:transparent,
                background_color_inside=:transparent, background_color_outside=:transparent)
            p3 = Plots.plot(t, power, xlabel="Time (s)", ylabel="Power (W)",
                label="Power", background_color=:transparent,
                background_color_inside=:transparent, background_color_outside=:transparent)
            p = Plots.plot(p1, p2, p3, layout=(3, 1),
                background_color=:transparent, background_color_inside=:transparent,
                background_color_outside=:transparent)
            Plots.savefig(p, joinpath(figdir, "pole_motion.pdf"))
            display(p)
        catch err
            @warn "Plots.jl not available; skipping plot." err
        end
    end
end

@testset "pole motion gradients" begin
    Rhub = 0.1
    Rtip = 1.0
    B = 2
    precone = 0.0
    rotor = CCBlade.Rotor(Rhub, Rtip, B, precone=precone, turbine=true)

    r = [0.2, 0.6, 0.9]
    chord = fill(0.1, 3)
    theta = fill(0.0, 3)

    function affun(alpha, Re, M)
        return 6.2 * alpha, 0.01
    end

    sections = [CCBlade.Section(r[i], chord[i], theta[i], affun) for i in eachindex(r)]

    Vhub = 10.0
    Omega = 20.0
    pitch = 0.0
    yaw = 0.0
    tilt = 0.0
    azimuth = 0.0
    hubHt = 10.0
    shearExp = 0.0
    rho = 1.0

    nsteps = 20
    dt = 0.1
    omega_m = 2.0 * pi / (nsteps * dt)

    function pole_motion_wrapper(x)
        arm_z = x[1]
        heave_amp = x[2]
        zero_val = zero(arm_z) + zero(heave_amp)

        roll_amp = 1.0 * pi / 180 + zero_val
        pitch_amp = 2.0 * pi / 180 + zero_val
        yaw_amp = 0.5 * pi / 180 + zero_val

        arm = (zero_val, zero_val, arm_z)
        params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.3, 3.0)

        heave0 = heave_amp * sin(zero_val)
        heave_dot0 = heave_amp * omega_m * cos(zero_val)
        base_pos0 = (zero_val, zero_val, heave0)
        base_vel0 = (zero_val, zero_val, heave_dot0)
        base_angles0 = (roll_amp * sin(zero_val), pitch_amp * sin(zero_val),
            yaw_amp * sin(zero_val))
        base_omega0 = (roll_amp * omega_m * cos(zero_val),
            pitch_amp * omega_m * cos(zero_val),
            yaw_amp * omega_m * cos(zero_val))

        ops0, info0 = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
            yaw, tilt, azimuth, hubHt, shearExp, rho;
            base_pos=base_pos0,
            base_vel=base_vel0,
            base_angles=base_angles0,
            base_omega=base_omega0,
            arm=arm)
        state = UnsteadyKineticRotorDynamics.UnsteadyState(sections, ops0; V_wake_old=info0.Vhub_eff)

        T_sum = zero_val
        for i in 1:nsteps
            t = (i - 1) * dt
            heave = heave_amp * sin(omega_m * t)
            heave_dot = heave_amp * omega_m * cos(omega_m * t)

            roll = roll_amp * sin(omega_m * t)
            pitch_b = pitch_amp * sin(omega_m * t)
            yaw_b = yaw_amp * sin(omega_m * t)

            roll_dot = roll_amp * omega_m * cos(omega_m * t)
            pitch_dot = pitch_amp * omega_m * cos(omega_m * t)
            yaw_dot = yaw_amp * omega_m * cos(omega_m * t)

            base_pos = (zero_val, zero_val, heave)
            base_vel = (zero_val, zero_val, heave_dot)
            base_angles = (roll, pitch_b, yaw_b)
            base_omega = (roll_dot, pitch_dot, yaw_dot)

            ops, _ = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
                yaw, tilt, azimuth, hubHt, shearExp, rho;
                base_pos=base_pos,
                base_vel=base_vel,
                base_angles=base_angles,
                base_omega=base_omega,
                arm=arm)

            out = UnsteadyKineticRotorDynamics.unsteady_step!(state, rotor, sections, ops, params; dt=dt)
            T, _ = CCBlade.thrusttorque(rotor, sections, out)
            T_sum += T
        end

        return T_sum / nsteps
    end

    x = [hubHt, 0.2]
    g = ForwardDiff.gradient(pole_motion_wrapper, x)
    g_fd = FiniteDiff.finite_difference_gradient(pole_motion_wrapper, x, Val{:central})

    @test all(isfinite.(g))
    @test all(isfinite.(g_fd))
    @test isapprox(g[1], g_fd[1]; rtol=5e-2, atol=5e-2)
    @test isapprox(g[2], g_fd[2]; rtol=5e-2, atol=5e-2)
end

@testset "pole motion validation" begin
    Rhub = 0.1
    Rtip = 1.0
    B = 2
    precone = 0.0

    r = [0.2, 0.6, 0.9]

    Vhub = 10.0
    Omega = 20.0
    pitch = 0.0
    yaw = 0.0
    tilt = 0.0
    azimuth = 0.0
    hubHt = 10.0
    shearExp = 0.0
    rho = 1.0

    arm = (0.0, 0.0, hubHt)

    function op_components(ops)
        Vx = [ops[i].Vx for i in eachindex(ops)]
        Vy = [ops[i].Vy for i in eachindex(ops)]
        return Vx, Vy
    end

    # zero-motion invariance
    ops_zero, info_zero = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        yaw, tilt, azimuth, hubHt, shearExp, rho; arm=arm)
    ops_ref = [CCBlade.windturbine_op(Vhub, Omega, pitch, r[i], precone,
        yaw, tilt, azimuth, hubHt, shearExp, rho) for i in eachindex(r)]

    Vx_zero, Vy_zero = op_components(ops_zero)
    Vx_ref, Vy_ref = op_components(ops_ref)
    @test all(isapprox.(Vx_zero, Vx_ref; atol=1e-12))
    @test all(isapprox.(Vy_zero, Vy_ref; atol=1e-12))
    @test isapprox(info_zero.Vhub_eff, Vhub; atol=1e-12)
    @test isapprox(info_zero.hubHt_eff, hubHt; atol=1e-12)
    @test isapprox(info_zero.yaw_eff, yaw; atol=1e-12)
    @test isapprox(info_zero.tilt_eff, tilt; atol=1e-12)

    # pure heave (position + velocity)
    heave = 0.5
    heave_dot = -0.2
    ops_heave, info_heave = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        yaw, tilt, azimuth, hubHt, shearExp, rho;
        base_pos=(0.0, 0.0, heave),
        base_vel=(0.0, 0.0, heave_dot),
        arm=arm)
    Vrel_x = Vhub
    Vrel_y = 0.0
    Vrel_z = -heave_dot
    Vhub_eff = sqrt(Vrel_x * Vrel_x + Vrel_y * Vrel_y + Vrel_z * Vrel_z)
    wind_tilt = atan(Vrel_z, hypot(Vrel_x, Vrel_y))
    tilt_eff = tilt - wind_tilt
    hubHt_eff = hubHt + heave
    ops_heave_ref = [CCBlade.windturbine_op(Vhub_eff, Omega, pitch, r[i], precone,
        yaw, tilt_eff, azimuth, hubHt_eff, shearExp, rho) for i in eachindex(r)]
    Vx_heave, Vy_heave = op_components(ops_heave)
    Vx_heave_ref, Vy_heave_ref = op_components(ops_heave_ref)
    @test all(isapprox.(Vx_heave, Vx_heave_ref; rtol=1e-10, atol=1e-10))
    @test all(isapprox.(Vy_heave, Vy_heave_ref; rtol=1e-10, atol=1e-10))
    @test isapprox(info_heave.hubHt_eff, hubHt_eff; atol=1e-12)
    @test isapprox(info_heave.Vhub_eff, Vhub_eff; rtol=1e-12, atol=1e-12)
    @test isapprox(info_heave.tilt_eff, tilt_eff; rtol=1e-12, atol=1e-12)

    # pure heave position only should not change inflow magnitude or yaw/tilt
    _, info_heave_pos = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r,
        precone, yaw, tilt, azimuth, hubHt, shearExp, rho;
        base_pos=(0.0, 0.0, heave),
        base_vel=(0.0, 0.0, 0.0),
        arm=arm)
    @test isapprox(info_heave_pos.Vhub_eff, Vhub; atol=1e-12)
    @test isapprox(info_heave_pos.yaw_eff, yaw; atol=1e-12)
    @test isapprox(info_heave_pos.tilt_eff, tilt; atol=1e-12)
    @test isapprox(info_heave_pos.hubHt_eff, hubHt + heave; atol=1e-12)

    # constant sway velocity updates yaw and inflow magnitude
    sway = 0.4
    _, info_sway = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        yaw, tilt, azimuth, hubHt, shearExp, rho;
        base_vel=(0.0, sway, 0.0),
        arm=arm)
    Vhub_eff_sway = sqrt(Vhub * Vhub + sway * sway)
    yaw_eff_sway = yaw + atan(sway, Vhub)
    @test isapprox(info_sway.Vhub_eff, Vhub_eff_sway; rtol=1e-12, atol=1e-12)
    @test isapprox(info_sway.yaw_eff, yaw_eff_sway; rtol=1e-12, atol=1e-12)

    # rigid-body rotation about base (omega about y, arm along z)
    omega_y = 0.2
    _, info_rot = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        yaw, tilt, azimuth, hubHt, shearExp, rho;
        base_omega=(0.0, omega_y, 0.0),
        arm=arm)
    hub_vel_expected = (omega_y * hubHt, 0.0, 0.0)
    @test isapprox(info_rot.hub_vel[1], hub_vel_expected[1]; atol=1e-12)
    @test isapprox(info_rot.hub_vel[2], hub_vel_expected[2]; atol=1e-12)
    @test isapprox(info_rot.hub_vel[3], hub_vel_expected[3]; atol=1e-12)
    @test isapprox(info_rot.Vhub_eff, abs(Vhub - hub_vel_expected[1]); rtol=1e-12, atol=1e-12)

    # small-angle linearization (base angles)
    eps = 1e-4
    _, info_small = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        yaw, tilt, azimuth, hubHt, shearExp, rho;
        base_angles=(eps, -2 * eps, 3 * eps),
        arm=arm)
    @test isapprox(info_small.yaw_eff, yaw + 3 * eps; atol=1e-8)
    @test isapprox(info_small.tilt_eff, tilt - 2 * eps; atol=1e-8)
    @test isapprox(info_small.azimuth_eff, azimuth + eps; atol=1e-8)
    @test isapprox(info_small.Vhub_eff, Vhub; atol=1e-10)

    # symmetry checks
    heave_dot_pos = 0.3
    _, info_sym_p = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        0.0, 0.0, azimuth, hubHt, shearExp, rho;
        base_vel=(0.0, 0.0, heave_dot_pos),
        arm=arm)
    _, info_sym_n = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        0.0, 0.0, azimuth, hubHt, shearExp, rho;
        base_vel=(0.0, 0.0, -heave_dot_pos),
        arm=arm)
    @test isapprox(info_sym_p.Vhub_eff, info_sym_n.Vhub_eff; atol=1e-12)
    @test isapprox(info_sym_p.tilt_eff, -info_sym_n.tilt_eff; atol=1e-10)

    sway_mag = 0.25
    _, info_sway_p = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        0.0, 0.0, azimuth, hubHt, shearExp, rho;
        base_vel=(0.0, sway_mag, 0.0),
        arm=arm)
    _, info_sway_n = UnsteadyKineticRotorDynamics.windturbine_op_motion(Vhub, Omega, pitch, r, precone,
        0.0, 0.0, azimuth, hubHt, shearExp, rho;
        base_vel=(0.0, -sway_mag, 0.0),
        arm=arm)
    @test isapprox(info_sway_p.yaw_eff, -info_sway_n.yaw_eff; atol=1e-10)
end
