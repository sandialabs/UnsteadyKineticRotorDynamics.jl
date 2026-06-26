using Test

import ForwardDiff
import FiniteDiff
import CCBlade
import UnsteadyKineticRotorDynamics

# ForwardDiff vs FiniteDiff check for time-averaged response.

@testset "windturbine time gradients" begin
    data_dir = joinpath(@__DIR__, "..", "data")

    # --- rotor definition ---
    Rhub = 1.5
    Rtip = 63.0
    B = 3
    precone = 2.5 * pi / 180

    # --- section definitions ---
    r = [2.8667, 5.6000, 8.3333, 11.7500, 15.8500, 19.9500, 24.0500,
        28.1500, 32.2500, 36.3500, 40.4500, 44.5500, 48.6500, 52.7500,
        56.1667, 58.9000, 61.6333]
    chord = [3.542, 3.854, 4.167, 4.557, 4.652, 4.458, 4.249, 4.007, 3.748,
        3.502, 3.256, 3.010, 2.764, 2.518, 2.313, 2.086, 1.419]
    theta = pi / 180 * [13.308, 13.308, 13.308, 13.308, 11.480, 10.162, 9.011, 7.795,
        6.544, 5.361, 4.188, 3.125, 2.319, 1.526, 0.863, 0.370, 0.106]

    # --- airfoils ---
    aftypes = Array{CCBlade.AlphaAF}(undef, 8)
    aftypes[1] = CCBlade.AlphaAF(joinpath(data_dir, "Cylinder1.dat"), radians = false)
    aftypes[2] = CCBlade.AlphaAF(joinpath(data_dir, "Cylinder2.dat"), radians = false)
    aftypes[3] = CCBlade.AlphaAF(joinpath(data_dir, "DU40_A17.dat"), radians = false)
    aftypes[4] = CCBlade.AlphaAF(joinpath(data_dir, "DU35_A17.dat"), radians = false)
    aftypes[5] = CCBlade.AlphaAF(joinpath(data_dir, "DU30_A17.dat"), radians = false)
    aftypes[6] = CCBlade.AlphaAF(joinpath(data_dir, "DU25_A17.dat"), radians = false)
    aftypes[7] = CCBlade.AlphaAF(joinpath(data_dir, "DU21_A17.dat"), radians = false)
    aftypes[8] = CCBlade.AlphaAF(joinpath(data_dir, "NACA64_A17.dat"), radians = false)

    # indices correspond to which airfoil is used at which station
    af_idx = [1, 1, 2, 3, 4, 4, 5, 6, 6, 7, 7, 8, 8, 8, 8, 8, 8]

    # create airfoil array
    airfoils = aftypes[af_idx]

    # define sections
    sections = [CCBlade.Section(r[i], chord[i], theta[i], airfoils[i])
                for i in eachindex(r)]

    # --- operating conditions ---
    yaw = 0.0 * pi / 180
    tilt = 5.0 * pi / 180
    hubHt = 90.0
    shearExp = 0.2

    Vinf_nom = 10.0
    tsr = 7.55
    pitch = 0.0
    rho = 1.225

    azangles = pi / 180 * [0.0, 90.0, 180.0, 270.0]
    naz = length(azangles)

    # --- time stepping ---
    t_final = 20.0
    dt = 0.5
    t = 0:dt:t_final
    nt = length(t)

    Vinf = Vinf_nom .* (0.75 .+ 0.25 .* sin.(2 * pi * t ./ t_final))

    # parameters that passthrough: airfoils, B
    function windturbinewrapper(x)
        # unpack
        Rtip = x[1]

        Rhubp, Rtipp = promote(Rhub, Rtip)
        preconep = convert(typeof(Rtipp), precone)
        rotorp = CCBlade.Rotor(Rhubp, Rtipp, B, precone = preconep, turbine = true)
        rotorR = Rtipp * cos(preconep)
        Omega = Vinf_nom * tsr / rotorR
        params = UnsteadyKineticRotorDynamics.UnsteadyParams(0.5, 2.0)
        function op_vec(V, az)
            [CCBlade.windturbine_op(V, Omega, pitch, r[i], preconep, yaw, tilt,
                 az, hubHt, shearExp, rho) for i in eachindex(r)]
        end

        zero_val = zero(Rtip)
        zero_out = CCBlade.Outputs(ntuple(_ -> zero_val, 15)...)

        states = Vector{UnsteadyKineticRotorDynamics.UnsteadyState}(undef, naz)
        for j in 1:naz
            op0 = op_vec(Vinf[1], azangles[j])
            states[j] = UnsteadyKineticRotorDynamics.UnsteadyState(sections, op0;
                init_output = zero_out,
                V_wake_old = Vinf[1] + zero_val)
        end

        cp_sum = zero_val
        ct_sum = zero_val

        for i in 1:nt
            V = Vinf[i]
            Tu = zero_val
            Qu = zero_val

            for j in 1:naz
                ops = op_vec(V, azangles[j])
                out = UnsteadyKineticRotorDynamics.unsteady_step!(
                    states[j], rotorp, sections, ops, params; dt = dt)
                Tj, Qj = CCBlade.thrusttorque(rotorp, sections, out)
                Tu += Tj / naz
                Qu += Qj / naz
            end

            cp, ct, _ = CCBlade.nondim(Tu, Qu, V, Omega, rho, rotorp, "windturbine")
            cp_sum += cp
            ct_sum += ct
        end

        return [cp_sum / nt; ct_sum / nt]
    end

    x = [Rtip]

    J = ForwardDiff.jacobian(windturbinewrapper, x)
    J_fd = FiniteDiff.finite_difference_jacobian(windturbinewrapper, x, Val{:central})

    @test all(isfinite.(J))
    @test all(isfinite.(J_fd))
    @test isapprox(J[1, 1], J_fd[1, 1]; rtol = 1e-2, atol = 1e-2)
    @test isapprox(J[2, 1], J_fd[2, 1]; rtol = 1e-2, atol = 1e-2)
end
