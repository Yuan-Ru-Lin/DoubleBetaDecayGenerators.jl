module DoubleBetaDecayGenerators

export ZeroNuDBDData, TwoNuDBDData, DBDGenerator

using DelimitedFiles
using StatsBase
using EmpiricalDistributions
using LinearAlgebra
using Rotations
using Random: AbstractRNG, default_rng
using Artifacts

"""
    acute_inv_cdf(x, q)

The inverse of the cumulative of the "acute" probability distribution, where q
is the linear angular coefficient.

[1] https://stats.stackexchange.com/questions/171592/generate-random-numbers-with-linear-distribution
"""
acute_inv_cdf(x::Real, q::Real) = (√(q^2 - 2q + 4q*x+1) - 1)/q

const DBD_ISOTOPES::NTuple{19,Symbol} = (
    :Ca48, :Cd116, :Gd160, :Ge76, :Mo100, :Nd148, :Nd150, :Pd110, :Pt198,
    :Se82, :Sm154, :Sn124, :Te128, :Te130, :Th232, :U238, :Xe134, :Xe136, :Zr96
)
abstract type DBDData end
abstract type DBDDataSource end
abstract type KotilaIachello2012 <: DBDDataSource end

"""
    ZeroNuDBDData(path)

A structure to hold numerical calculations of zero-neutrino double-beta decay in
the format provided by Iachello, F. and Kotila, J. [1].

[1] https://nucleartheory.yale.edu/double-beta-decay-phase-space-factors
"""
struct ZeroNuDBDData <: DBDData
    _cor_data::AbstractArray{Real, 1}
    ses_dist::UvBinnedDist
    Q_value_keV::Real
    function ZeroNuDBDData(isotope::Symbol, datasource::Type{<:DBDDataSource} = KotilaIachello2012)
        @assert isotope in DBD_ISOTOPES "$isotope does not undergo double-beta decay.\nAvailable isotopes are $DBD_ISOTOPES."
        @debug "Reading raw data from '$path'"
        path = @artifact_str("$(nameof(datasource))/$isotope")

        _cor_data = readdlm("$path/cor_0v.txt", Float64)[:, 3]
        _ses_data = readdlm("$path/ses_0v.txt", Float64)[:, 3]
        e_max_keV = length(_ses_data)

        ses_dist = UvBinnedDist(Histogram(0:e_max_keV, _ses_data, :left, true))

        new(_cor_data, ses_dist, Float64(e_max_keV))
    end
end

function Base.rand(rng::AbstractRNG, data::ZeroNuDBDData)
    # sample electron energies
    E1 = rand(rng, data.ses_dist)
    E2 = data.Q_value_keV - E1
    @debug "Electron energies" E1, E2

    # sample the angle between the two electrons
    # select angular coefficient
    idx = Int(trunc(E1)) # the energy in keV can be used as index
    α = data._cor_data[idx > 0 ? idx : 1]
    # sample cosine of the angle from a linear pdf
    cosθ12 = acute_inv_cdf(rand(rng), α)

    E1, E2, cosθ12
end

"""
    TwoNuDBDData(path)

A structure to hold numerical calculations of two-neutrino double-beta decay in
the format provided by Iachello, F. and Kotila, J. [1].

[1] https://nucleartheory.yale.edu/double-beta-decay-phase-space-factors
"""
struct TwoNuDBDData <: DBDData
    _ses_data::AbstractArray{Real, 1}
    _cor_data::AbstractArray{Real, 1}
    _tds_data::AbstractArray{Real, 2}
    ses_dist::UvBinnedDist
    tds_dist::MvBinnedDist

    function TwoNuDBDData(isotope::Symbol, datasource::Type{<:DBDDataSource} = KotilaIachello2012)
        @assert isotope in DBD_ISOTOPES "$isotope does not undergo double-beta decay.\nAvailable isotopes are $DBD_ISOTOPES."
        @debug "Reading raw data from '$path'"
        path = @artifact_str("$(nameof(datasource))/$isotope")

        _ses_data = readdlm("$path/ses.txt", Float64)[:, 3]
        e_max_keV = length(_ses_data)
        _cor_data = readdlm("$path/cor.txt", Float64)[:, 3]
        raw_data  = readdlm("$path/2ds.txt", Float64)

        _tds_data = zeros(e_max_keV, e_max_keV)
        for i in 1:size(raw_data, 1)
            _tds_data[Int(raw_data[i, 1]), Int(raw_data[i, 2])] = raw_data[i, 5]
            _tds_data[Int(raw_data[i, 2]), Int(raw_data[i, 1])] = raw_data[i, 5]
        end

        ses_dist = UvBinnedDist(Histogram(0:e_max_keV, _ses_data, :left, true))
        tds_dist = MvBinnedDist(Histogram((0:e_max_keV, 0:e_max_keV), _tds_data, :left, true))

        new(_ses_data, _cor_data, _tds_data, ses_dist, tds_dist)
    end
end

function Base.rand(rng::AbstractRNG, data::TwoNuDBDData)
    # sample electron energies
    E1, E2 = rand(rng, data.tds_dist)
    @debug "Electron energies" E1, E2

    # sample the angle between the two electrons
    # select angular coefficient
    idx = Int(trunc(E1)) # the energy in keV can be used as index
    α = data._cor_data[idx > 0 ? idx : 1]
    # sample cosine of the angle from a linear pdf
    cosθ12 = acute_inv_cdf(rand(rng), α)
    @debug "Angle between the two electrons" cosθ12

    E1, E2, cosθ12
end

Base.rand(data::T) where {T<:DBDData} = rand(default_rng(), data)

"""
Legacy code for generating events in the DECAY0 format.
"""
module Legacy

using ..DoubleBetaDecayGenerators: TwoNuDBDData, rand, acute_inv_cdf
using Rotations
using ProgressMeter
using Printf

"""
Convert energy to momentum (relativistic)
"""
momentum(E) = √((E + 511)^2 - 511^2)

"""
    rand_2nbb(data::TwoNuDBDData)

Generate the momenta of the two electrons emitted in two-neutrino double-beta
decay, according to numerical calculations wrapped by 'data'. Units of keV.
"""
function rand_2nbb(data::TwoNuDBDData)

    # sample electron energies
    E1, E2 = rand(data.tds_dist)
    @debug "Electron energies" E1, E2

    # sample the angle between the two electrons
    # select angular coefficient
    idx = Int(trunc(E1)) # the energy in keV can be used as index
    α = data._cor_data[idx > 0 ? idx : 1]
    # sample cosine of the angle from a linear pdf
    cosθ12 = acute_inv_cdf(rand(), α)
    ϕ = 2π * rand()
    @debug "Direction of second electron, when the first is along z" cosθ12, ϕ

    p1 = [0, 0, 1]
    p2 = [cos(ϕ), sin(ϕ), cosθ12] ./ √(1 + cosθ12^2)

    # rotate randomly the two vectors with the same transformation,
    # to preserve the angle between them
    R = rand(RotMatrix{3})
    p1 = R * p1
    p2 = R * p2

    # transform it to momentum
    p1 = p1 .* momentum(E1)
    p2 = p2 .* momentum(E2)

    return p1, p2
end

"""
    dk0gen(N; options)

Generate N events and dump them on file in the DECAY0 format.
"""
function dk0gen(n::Int64; output="ge76-ssd.dk0", input="76Ge_ssd")

    data = TwoNuDBDData(input)

    @info "Generating and writing events to file '$output'"
    fout = open(output, "w")

    print(fout,"""
 DECAY0 generated file: $output

  date and time :   12. 9.2020    17:50:37
  initial scrolling of (0,1) random number generator :            0

  event type: Ge76
              2nubb     0+ -> 0+     {2n}
              level, Elevel (MeV) =  0+       0.0000     MeV

 Format of data:
  for each event    - event's number,
                      time of event's start,
                      number of emitted particles;
  for each particle - GEANT number of particle,
                      x,y,z components of momentum,
                      time shift from previous time

 Time - in sec, momentum - in MeV/c

 First event and full number of events:
           1    $n

""")

    @showprogress for i in 1:n
        println(fout, "       $i  0.00000       2")

        p1, p2 = rand_2nbb(data)
        p1 = p1 ./ 1e3
        p2 = p2 ./ 1e3

        @printf fout "  3 %9.6f     %9.6f     %9.6f     0.00000\n" p1[1] p1[2] p1[3]
        @printf fout "  3 %9.6f     %9.6f     %9.6f     0.00000\n" p2[1] p2[2] p2[3]
    end

    close(fout)
end

end # module Legacy

function DBDGenerator end

end # module
