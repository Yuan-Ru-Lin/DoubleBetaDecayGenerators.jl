using Test
using DoubleBetaDecayGenerators
using Artifacts
using Random

@testset "Test control of random seed" begin
    rng = Xoshiro()

    Random.seed!(rng, 999)

    dat_2nu = TwoNuDBDData(:Ge76)
    E1, E2, cosθ12 = rand(rng, dat_2nu)
    @test E1 ≈ 374.2355832686694 atol=1e-15
    @test E2 ≈ 1251.059139039826 atol=1e-15
    @test cosθ12 ≈ 0.9180839309783714 atol=1e-15

    Random.seed!(rng, 999)

    dat_0nu = ZeroNuDBDData(:Ge76)
    E1, E2, cosθ12 = rand(rng, dat_0nu)
    @test E1 ≈ 2016.2355832686694 atol=1e-15
    @test E2 ≈ 22.764416731330584 atol=1e-15
    @test cosθ12 ≈ -0.9113939723175647 atol=1e-15
end
