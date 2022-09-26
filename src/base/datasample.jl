abstract type DataSample  end

struct Spectrum <: DataSample
    data::Vector
end

struct RealSpectrum <: DataSample 
    data::Vector{T} where T<:Real
end

struct UnfoldedSpectrum <: DataSample 
    data::Vector{T} where T<:Real
end

struct ComplexSpectrum <: DataSample 
    data::Vector{T} where T<:Complex
end

struct SpectralEnsemble 
    data::Vector{T} where T<:DataSample
end 