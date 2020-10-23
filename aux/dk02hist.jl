using StatsBase
using LinearAlgebra
using Plots
using CSV
using Tables

events = []

open("ge76-hsd.dk0") do f
# open("2nbb-decay0.txt") do f

    # skip header
    for _ in 1:22 readline(f) end

    while !eof(f)

        # skip first line
        readline(f)

        # register new event
        push!(events, [])

        # parse momentum
        p = split(readline(f), " ", keepempty=false)
        push!(last(events), (parse(Float64, p[2]), parse(Float64, p[3]), parse(Float64, p[4])))
        p = split(readline(f), " ", keepempty=false)
        push!(last(events), (parse(Float64, p[2]), parse(Float64, p[3]), parse(Float64, p[4])))

    end
end

function kine(p)
    return 1E3(sqrt(p[1]^2 + p[2]^2 + p[3]^2 + 0.511^2) - 0.511)
end

function costheta(p1, p2)
    dot(p1, p2)/norm(p1)/norm(p2)
end

sums = Number[]
corr = Number[]
ses = Number[]

for e in events
    push!(sums, sum(kine, e))
    push!(corr, costheta(e[1], e[2]))
    push!(ses, kine(e[1]))
end

h_hsd = normalize(fit(Histogram, sums, 0:10:2039, closed=:left))
plot!(h_hsd, st=:step)

# h = normalize(fit(Histogram, corr[(ses .< 300) .& (ses .> 200)], nbins=20))
# plot(h, st=:step, xlabel="cos(theta)")

# h2 = fit(Histogram, (ses, corr))
# plot(h2, xlabel="energy E_1", ylabel="cos(theta_12)")

# plot(sums[1:10000], corr[1:10000], st=:scatter)

# CSV.write("dk0_save.csv", Tables.table([collect(h.edges)[1][1:end-1] h.weights]))
