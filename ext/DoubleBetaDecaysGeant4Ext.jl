module DoubleBetaDecaysGeant4Ext

using DoubleBetaDecays
using Geant4
using Geant4.SystemOfUnits
using StaticArrays
using Rotations

@kwdef mutable struct ZeroNuDBDGeneratorData <: G4JLGeneratorData
    gun::Union{Nothing,CxxPtr{G4ParticleGun}} = nothing
    spectrum::ZeroNuDBDData = ZeroNuDBDData()
end

function DoubleBetaDecays.ZeroNuDBDGenerator(; kwargs...)

    data = ZeroNuDBDGeneratorData(; kwargs...)

    function _init(data::ZeroNuDBDGeneratorData, app::G4JLApplication)::Nothing
        data.gun = move!(G4ParticleGun())
        # XXX: Why doesn't this work here?
        #SetParticleDefinition(data.gun, FindParticle("electron"))
        return nothing
    end

    function _gen(evt::G4Event, data::ZeroNuDBDGeneratorData)::Nothing
        # XXX: Where do I store `rng` if I want to control the random seed?
        E1, E2, cosθ12 = rand(data.spectrum)

        sinθ12 = sqrt(1.0 - cosθ12^2)
        ϕ = 2π * rand()
        dir1 = SVector{3}([0.0, 0.0, 1.0])
        dir2 = SVector{3}(sinθ12 * cos(ϕ), sinθ12 * sin(ϕ), cosθ12)
        rot = rand(RotMatrix{3})
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

    G4JLPrimaryGenerator("ZeroNuDBDGenerator", data; init_method=_init, generate_method=_gen)
end

@kwdef mutable struct TwoNuDBDGeneratorData <: G4JLGeneratorData
    gun::Union{Nothing,CxxPtr{G4ParticleGun}} = nothing
    spectrum::TwoNuDBDData = TwoNuDBDData()
end

function DoubleBetaDecays.TwoNuDBDGenerator(; kwargs...)

    data = TwoNuDBDGeneratorData(; kwargs...)

    function _init(data::TwoNuDBDGeneratorData, app::G4JLApplication)::Nothing
        data.gun = move!(G4ParticleGun())
        # XXX: Why doesn't this work here?
        #SetParticleDefinition(data.gun, FindParticle("electron"))
        return nothing
    end

    function _gen(evt::G4Event, data::TwoNuDBDGeneratorData)::Nothing
        # XXX: Where do I store `rng` if I want to control the random seed?
        E1, E2, cosθ12 = rand(data.spectrum)

        sinθ12 = sqrt(1.0 - cosθ12^2)
        ϕ = 2π * rand()
        dir1 = SVector{3}([0.0, 0.0, 1.0])
        dir2 = SVector{3}(sinθ12 * cos(ϕ), sinθ12 * sin(ϕ), cosθ12)
        rot = rand(RotMatrix{3})
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

    G4JLPrimaryGenerator("TwoNuDBDGenerator", data; init_method=_init, generate_method=_gen)
end

end
