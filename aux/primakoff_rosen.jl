using Random
using Distributions

struct DoubleBetaSumEnergyPRContinuous <: ContinuousUnivariateDistribution
    Q_ββ::Float64
    function DoubleBetaSumEnergyPRContinuous(Q_ββ::T) where {T<:Real}
        Distributions.@check_args(DoubleBetaSumEnergyPRContinuous, Q_ββ >= 0)
        new(Q_ββ)
    end
end

Distributions.params(d::DoubleBetaSumEnergyPRContinuous) = d.Q_ββ
Distributions.minimum(::DoubleBetaSumEnergyPRContinuous) = 0
Distributions.maximum(d::DoubleBetaSumEnergyPRContinuous) = d.Q_ββ

function Distributions.pdf(d::DoubleBetaSumEnergyPRContinuous, x::T) where {T<:Real}
    m = 501.998
    x /= m
    Q = d.Q_ββ / m

    expr = x * (x^4 + 10x^3 + 40x^2 + 60x + 30) * (Q - x)^5
    expr /= Q^7 * (1980 + 990Q + 220Q^2 + 22Q^3 + Q^4) / 2772

    return expr
end

function Distributions.cdf(d::DoubleBetaSumEnergyPRContinuous, x::T) where {T<:Real}
    m = 501.998
    x /= m
    Q = d.Q_ββ / m

    expr = 15Q^5*x^2 + 10Q^4*(2Q-5)*x^3 + 5Q^3*(Q*(2Q-15) + 15)*x^4
    expr += 2Q^2*(Q*((Q-20)*Q + 60) - 30)*x^5 + 1/2*(Q-2)*x^10 - 10/9*(Q-4)*(Q-1)*x^9
    expr += 5/4*(Q*((Q-10)*Q + 20) - 6)*x^8 - 5/7*(Q*(Q*((Q-20)*Q + 80) - 60) + 6)*x^7
    expr += 1/6*Q*(Q*((Q - 40)*(Q - 10)*Q - 600) + 150)*x^6 - x^11/11

    expr /= Q^7 * (1980 + 990Q + 220Q^2 + 22Q^3 + Q^4) / 2772

    return expr
end

struct DoubleBetaSingleEnergyPRContinuous <: ContinuousUnivariateDistribution
    Q_ββ::Float64
    function DoubleBetaSingleEnergyPRContinuous(Q_ββ::T) where {T<:Real}
        Distributions.@check_args(DoubleBetaSingleEnergyPRContinuous, Q_ββ >= 0)
        new(Q_ββ)
    end
end

Distributions.params(d::DoubleBetaSingleEnergyPRContinuous) = d.Q_ββ
Distributions.minimum(::DoubleBetaSingleEnergyPRContinuous) = 0
Distributions.maximum(d::DoubleBetaSingleEnergyPRContinuous) = d.Q_ββ

function Distributions.pdf(d::DoubleBetaSingleEnergyPRContinuous, x::T) where {T<:Real}
    m = 501.998
    x /= m
    Q = d.Q_ββ / m

    expr = (x+1)^2 * (Q-x)^6 * ((Q-x)^2 + 8*(Q-x) + 28)
    expr /= (Q^7 * (1980 + 990Q + 220Q^2 + 22Q^3 + Q^4))/495

    return expr
end

struct DoubleBetaPR
    Q_ββ::Float64
    sum_energy::DiscreteNonParametric
    single_energy::DiscreteNonParametric

    function DoubleBetaPR(Q_ββ::T) where {T<:Real}

        nsamples = 1E5
        xs_sums = Float64[]
        xp_sums = Float64[]
        sums_pr = DoubleBetaSumEnergyPRContinuous(Q_ββ)

        for i in 1:nsamples
            push!(xs_sums, (i-1)*(Distributions.maximum(sums_pr)-Distributions.minimum(sums_pr))/nsamples)
            push!(xp_sums, Distributions.pdf(sums_pr, last(xs_sums)))
        end

        xp_sums = xp_sums ./ sum(xp_sums)

        xs_ses = Float64[]
        xp_ses = Float64[]
        ses_pr = DoubleBetaSingleEnergyPRContinuous(Q_ββ)

        for i in 1:nsamples
            push!(xs_ses, (i-1)*(Distributions.maximum(ses_pr)-Distributions.minimum(ses_pr))/nsamples)
            push!(xp_ses, Distributions.pdf(ses_pr, last(xs_ses)))
        end

        xp_ses = xp_ses ./ sum(xp_ses)

        new(Q_ββ,
            DiscreteNonParametric(xs_sums, xp_sums),
            DiscreteNonParametric(xs_ses, xp_ses))
    end
end
