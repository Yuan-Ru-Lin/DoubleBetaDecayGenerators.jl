using DoubleBetaDecayGenerators

function @main(args::Vector{String})::Cint
    dat_2nu = TwoNuDBDData()
    for i in 1:10
        println(Core.stdout, rand(dat_2nu))
    end
    return 0
end
