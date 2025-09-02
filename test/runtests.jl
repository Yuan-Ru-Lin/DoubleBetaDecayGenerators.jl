using Test
using DoubleBetaDecayGenerators
using Random

@testset "Test control of random seed" begin
    rng = Xoshiro()

    Random.seed!(rng, 999)

    dat_2nu = TwoNuDBDData()
    E1, E2, cosθ12 = rand(rng, dat_2nu)
    @test E1 ≈ 374.2355832686694 atol=1e-15
    @test E2 ≈ 1251.059139039826 atol=1e-15
    @test cosθ12 ≈ 0.9180839309783714 atol=1e-15

    # XXX: For ZeroNuDBDData the test fails with Julia v1.11.6
    #Random.seed!(rng, 999)

    #dat_0nu = ZeroNuDBDData()
    #E1, E2, cosθ12 = rand(rng, dat_0nu)
    #@test E1 ≈ 1628.9956960725929 atol=1e-15
    #@test E2 ≈ 410.00430392740714 atol=1e-15
    #@test cosθ12 ≈ -0.7239804777193068 atol=1e-15
end
