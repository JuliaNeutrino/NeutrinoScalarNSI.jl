using NeutrinoScalarNSI, GLoBES, Test

# P(μ → e), P(μ → μ) of the C engine of SNSI-Decoherence-DUNE: config, sign of CP, E, P_μe, P_μμ
const C_REFERENCE = [
    parse.(Float64, split(l)) for
    l in eachline(joinpath(@__DIR__, "snsi_c_reference.txt")) if !startswith(l, "#")
]

const DUNE = (; L = 1300, ρ = 2.848)
const E = [0.5, 1.0, 2.5, 5.0]

@testset "parameters" begin
    m = NewPhysics(ScalarNSI(:eμ, :tautau))
    @test m.names == [:m_lightest, :η_eμ, :δ_eμ, :η_ττ]
    @test length(NewPhysics(ScalarNSI()).names) == 1 + 6 + 3
    q = OscParams(m, nufit_params("6.1"); eta_emu = 0.1, delta_emu = 1.0, mlightest = 0.01)
    @test (q.η_eμ, q.δ_eμ, q.m_lightest, q.η_ττ) == (0.1, 1.0, 0.01, 0.0)
    @test_throws ArgumentError ScalarNSI(:ex)
    @test_throws ArgumentError ScalarNSI(:ee; scale = :Δm21)
    @test_throws ArgumentError OscParams(m; θ12 = 0.6, ε_ee = 0.1)
    @test sprint(show, ScalarNSI(:eτ)) == "ScalarNSI(eτ)"
    @test sprint(show, ScalarNSI(:ee; scale = 0.05)) == "ScalarNSI(ee; scale = 0.05)"
end

@testset "probabilities" begin
    p = nufit_params("6.1")
    standard = [probability(p, :μ, :e; E, DUNE..., antineutrino = a) for a in (false, true)]
    m = NewPhysics(ScalarNSI(:eμ, :eτ))
    use_model!(m)

    # η = 0: the standard probabilities, whatever the lightest mass
    for mℓ in (0.0, 0.05)
        q = OscParams(m, p; m_lightest = mℓ)
        @test probability(q, :μ, :e; E, DUNE...) ≈ standard[1] rtol = 1e-10
        @test probability(q, :μ, :e; E, DUNE..., antineutrino = true) ≈ standard[2] rtol = 1e-10
    end

    # The C engine (config 2: scalar NSI only)
    base = OscParams(
        θ12 = asin(sqrt(0.3088)),
        θ13 = asin(sqrt(0.02249)),
        θ23 = asin(sqrt(0.470)),
        δ = deg2rad(207),
        Δm21 = 7.537e-5,
        Δm31 = 2.521e-3,
    )
    q = OscParams(m, base; η_eμ = 0.05, δ_eμ = 1.0)
    for (config, cp, Eν, pμe, pμμ) in C_REFERENCE
        anti = cp < 0
        @test probability(q, :μ, :e; E = [Eν], DUNE..., antineutrino = anti)[1] ≈ pμe rtol = 1e-10
        @test probability(q, :μ, :μ; E = [Eν], DUNE..., antineutrino = anti)[1] ≈ pμμ rtol = 1e-10
    end

    # With η ≠ 0 the lightest mass matters
    q = OscParams(m, p; η_eτ = 0.1)
    @test !(
        probability(q, :μ, :e; E, DUNE...) ≈
        probability(OscParams(q; m_lightest = 0.05), :μ, :e; E, DUNE...)
    )

    # A fixed scale equal to √|Δm²31| is the default convention
    fixed = NewPhysics(ScalarNSI(:eμ, :eτ; scale = sqrt(abs(p.Δm31))))
    use_model!(fixed)
    @test probability(OscParams(fixed, q), :μ, :e; E, DUNE...) ≈
          (use_model!(m); probability(q, :μ, :e; E, DUNE...)) rtol = 1e-12

    # Probabilities stay unitary
    use_model!(m)
    q = OscParams(m, p; η_eμ = 0.1, δ_eμ = 0.7, η_eτ = -0.05, m_lightest = 0.02)
    total = sum(probability(q, :μ, β; E, DUNE...) for β in (:e, :μ, :τ))
    @test total ≈ ones(length(E)) rtol = 1e-12
    standard_engine!()
end
