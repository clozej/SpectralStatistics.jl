push!(LOAD_PATH,"../src/")
using Documenter, SpectralStatistics
makedocs(sitename="SpectralStatistics.jl",
pages = [
    "index.md",
    "tutorial.md",
    "API.md"
],
format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true"
)

)

deploydocs(
    repo = "github.com/clozej/SpectralStatistics.jl.git",
)
