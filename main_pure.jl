function (@main)(args)
    jobid = parse(Int, args[1])
    ζ = parse(Float64, args[2])
    D = parse(Int, args[3])

    χ = collect(16:25)[jobid]

    βc = 1 / 1.2737

    timestep = 0.02βc / 128

    method = TEBD(; truncalg=TimeEvolutionPEPO.SU(ξ=ζ))

    X, _, Z = PAULI

    model = Model(TimeEvolutionPEPO.Ising(-1.0, Z, [2.5 * LocalOp(X)]))

    sim = Simulation(
        model,
        fill(ThermalState(), 2, 2)
        ComplexSpace(2),
        ComplexSpace(D);
        timestep=timestep,
        method=method
    )


    time = Float64[]

    δ = Float64[]
    ξ = Float64[]
    m = Float64[]

    obsalg = VUMPS(; bonddim=χ, tol=1e-8, maxiter=200)

    try
        simulate!(identity, sim; numsteps=2400)
        simulate!(sim; numsteps=1600, maxshots=40, maxtime=20.0 * 60.0 * 60.0) do sim

            psi = quantumstate(sim)

            dm = PurifiedDensityMatrix(psi, obsalg)

            eps1, eps2 = TimeEvolutionPEPO._correlationlength(dm)

            rdmA = partialtrace(dm, (1, 1))
            rdmB = partialtrace(dm, (2, 1))

            mA = abs(expval(rdmA, Z))
            mB = abs(expval(rdmB, Z))

            push!(time, sim.info.iterations * 2 * timestep)
            push!(δ, abs(eps1 - eps2))
            push!(ξ, 2 / eps1)
            push!(m, 1 / 2 * (mA + mB))

            return nothing
        end
    catch _
        rethrow()
    finally

        psi = quantumstate(sim)
        dm = PurifiedDensityMatrix(psi, obsalg)

        jldsave("χ=$(χ)_ζ=$(ζ)_D=$(D)_obs.jld2"; time, psi, dm, sim.info, χ, ζ, ξ, δ, D, m)
        jldsave("χ=$(χ)_ζ=$(ζ)_D=$(D)_sim.jld2"; sim)
    end

    return 0
end
