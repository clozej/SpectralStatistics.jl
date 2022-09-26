include("../base/probability.jl")
include("../base/datasample.jl")

using StatsBase

function mode_fluctuations(spect::UnfoldedSpectrum)
    e = spect.data
    e0 = round(e[1]) - 1.0
    println(e0)
    s = [float(i) for i in 1:length(e)] 
    return e, s .- e .+ e0
end

#=
function mode_fluctuations(spect::UnfoldedSpectrum, pts)
    e = spect.data
    e0 = round(e[1])
    h = StatsBase.fit(Histogram, e, pts)
    return h.weights
end
=#