using Documenter, Literate, DoubleBetaDecayGenerators

Literate.markdown("src/examples/plots.jl", "src/examples", execute = true, documenter = true, credit = true)

makedocs(
    sitename = "DoubleBetaDecayGenerators.jl",
    modules = [DoubleBetaDecayGenerators, DoubleBetaDecayGenerators.Legacy],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://JuliaHEP.github.io/DoubleBetaDecayGenerators.jl/stable/",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
        #"Getting Started" => "getting_started.md",
        "Examples" => ["examples/plots.md",],
    ],
)

deploydocs(
    repo = "github.com/JuliaHEP/DoubleBetaDecayGenerators.jl.git",
    forcepush = true,
    push_preview = true,
)
