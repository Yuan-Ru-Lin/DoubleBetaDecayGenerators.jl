module DoubleBetaDecayGeneratorsGeant4Ext

using DoubleBetaDecayGenerators
using DoubleBetaDecayGenerators: DBDData
using Geant4
using Geant4.SystemOfUnits
using StaticArrays
using Rotations
using Random

@kwdef mutable struct DBDGeneratorData{T<:DBDData, U<:AbstractRNG} <: G4JLGeneratorData
    gun::Union{Nothing,CxxPtr{G4ParticleGun}} = nothing
    spectrum::T = T()
    rng::U = default_rng()
end

function DoubleBetaDecayGenerators.DBDGenerator(::Type{T}; kwargs...) where {T<:DBDData}

    data = DBDGeneratorData(; spectrum=T(), kwargs...)

    function _init(data::U, app::G4JLApplication)::Nothing where {U<:DBDGeneratorData}
        data.gun = move!(G4ParticleGun())
        SetParticleDefinition(data.gun, FindParticle("e-"))
        return nothing
    end

    function _gen(evt::G4Event, data::U)::Nothing where {U<:DBDGeneratorData}
        E1, E2, cosθ12 = rand(data.rng, data.spectrum)

        sinθ12 = sqrt(1.0 - cosθ12^2)
        ϕ = 2π * rand(data.rng)
        dir1 = SVector{3}([0.0, 0.0, 1.0])
        dir2 = SVector{3}(sinθ12 * cos(ϕ), sinθ12 * sin(ϕ), cosθ12)
        rot = rand(data.rng, RotMatrix{3})
        dir1 = rot * dir1
        dir2 = rot * dir2

        # TODO: Implement decay position sampling within the detector volume
        #SetParticlePosition(data.gun, CxxPtr(G4ThreeVector(0.0, 0.0, 0.0)))

        SetParticleEnergy(data.gun, E1 * keV)
        SetParticleMomentumDirection(data.gun, G4ThreeVector(dir1...))
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))

        SetParticleEnergy(data.gun, E2 * keV)
        SetParticleMomentumDirection(data.gun, G4ThreeVector(dir2...))
        GeneratePrimaryVertex(data.gun, CxxPtr(evt))

        return nothing
    end

    G4JLPrimaryGenerator(string(T), data; init_method=_init, generate_method=_gen)
end

end
