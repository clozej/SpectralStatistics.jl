module SpectralStatistics

using Polynomials, LsqFit 
using StatsBase
using SpecialFunctions 
using Base.MathConstants
using InteractiveUtils

export DataSample
"""
    DataSample

`DataSample` is an abstract supertype encompassing all concrete implementations of 
data structures used for processing spectra in the SpectralStatistics.jl library.

"""
abstract type DataSample  end

include("base/datasample.jl")

include("base/probability.jl")
include("base/unfolding.jl")

export Model
"""
    Model

`Model` is an abstract supertype encompassing all concrete implementations of 
the analytical spectral statistics models in the SpectralStatistics.jl library.

"""
abstract type Model end

include("models/poisson.jl")
include("models/goe.jl")
include("models/gue.jl")
include("models/gse.jl")
include("models/brody.jl")
include("models/berryrobnik.jl")
include("models/berryrobnikbrody.jl")
include("models/gamma.jl")

include("statistics/modefluctuations.jl")
include("statistics/levelspacings.jl")
include("statistics/numbervariance.jl")
include("statistics/rigidity.jl")
include("statistics/spectralformfactor.jl")
include("statistics/lengthspectrum.jl")


include("utils/utils.jl")

end
