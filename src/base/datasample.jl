export RealSpectrum, UnfoldedSpectrum

#=
struct Spectrum <: DataSample
    data::Vector
end
=#

"""
    RealSpectrum <: DataSample

`RealSpectrum` is a concrete type used as a container for spectra of real numbers.

"""
struct RealSpectrum <: DataSample 
    data::Vector{T} where T<:Real
end


"""
    UnfoldedSpectrum <: DataSample

`UnfoldedSpectrum` is a concrete type used as a container for spectra after unfolding.

"""
struct UnfoldedSpectrum <: DataSample 
    data::Vector{T} where T<:Real
end

#=
struct ComplexSpectrum <: DataSample 
    data::Vector{T} where T<:Complex
end
=#
#=
struct SpectralEnsemble 
    data::Vector{T} where T<:DataSample
end
=# 