using DelimitedFiles
using StatsBase
using EmpiricalDistributions
using LinearAlgebra
using Rotations
using ProgressMeter
using Printf

"""
    TwoNuData(path)

A structure to hold numerical calculations of two-neutrino double-beta decay in
the format provided by Iachello, F. and Kotila, J. [1]. Must be initialized
with the path to the text files.

[1] https://nucleartheory.yale.edu/double-beta-decay-phase-space-factors
"""
struct TwoNuData
    _ses_data::AbstractArray{Real, 1}
    _cor_data::AbstractArray{Real, 1}
    _tds_data::AbstractArray{Real, 2}
    ses_dist::UvBinnedDist
    tds_dist::MvBinnedDist

    function TwoNuData(path::AbstractString)
        @info "Reading raw data from '$path'"

        _ses_data = readdlm("$path/ses.txt", Float64)[:, 3]
        _cor_data = readdlm("$path/cor.txt", Float64)[:, 3]
        raw_data  = readdlm("$path/2ds.txt", Float64)

        _tds_data = zeros(2039, 2039)
        for i in 1:size(raw_data, 1)
            _tds_data[Int(raw_data[i, 1]), Int(raw_data[i, 2])] = raw_data[i, 5]
            _tds_data[Int(raw_data[i, 2]), Int(raw_data[i, 1])] = raw_data[i, 5]
        end

        ses_dist = UvBinnedDist(Histogram(0:2039, _ses_data, :left, true))
        tds_dist = MvBinnedDist(Histogram((0:2039, 0:2039), _tds_data, :left, true))

        new(_ses_data, _cor_data, _tds_data, ses_dist, tds_dist)
    end
end

"""
    acute_inv_cdf(x, q)

The inverse of the cumulative of the "acute" probability distribution, where q
is the linear angular coefficient.

[1] https://stats.stackexchange.com/questions/171592/generate-random-numbers-with-linear-distribution
"""
acute_inv_cdf(x::Real, q::Real) = (√(q^2 - 2q + 4q*x+1) - 1)/q

"""
Convert energy to momentum (relativistic)
"""
momentum(E) = √((E + 511)^2 - 511^2)

"""
    rand_2nbb(data::TwoNuData)

Generate the momenta of the two electrons emitted in two-neutrino double-beta
decay, according to numerical calculations wrapped by 'data'. Units of keV.
"""
function rand_2nbb(data::TwoNuData)

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

global data = false

"""
    dk0gen(N; options)

Generate N events and dump them on file in the DECAY0 format.
"""
function dk0gen(n::Int64, file_name="ge76-ssd.dk0")

    if data == false
        global data = TwoNuData("76Ge_ssd")
    end

    @info "Generating and writing events to file '$file_name'"
    fout = open(file_name, "w")

    print(fout,"""
 DECAY0 generated file: $file_name

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
