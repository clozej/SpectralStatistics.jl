include("datasample.jl")

using Polynomials
using LsqFit


function unfold_spectrum(spect::DataSample, f::Function ) #f is the unfolding function
    unfolded = @. f(spect.data)
    return UnfoldedSpectrum(unfolded)
end

function unfold_spectrum(spect::DataSample, n::Int ) #f is the unfolding function
    coeffs = fit_integrated_density(spect, n)
    unfolded = [evalpoly(x, coeffs) for x in spect.data]
    return UnfoldedSpectrum(unfolded)
end

function fit_integrated_density(spect::DataSample, f::Function, p0) 
    x = spect.data
    y = [float(i) for i in 1:length(x)]
    fit = LsqFit.curve_fit(f, x, y, p0)
    return fit.param
end

function fit_integrated_density(spect::DataSample, n::Int) 
    x = spect.data
    y = [float(i) for i in 1:length(x)]
    fit = Polynomials.fit(x,y,n)
    return fit.coeffs
end

function weyl_law(spect::DataSample, A, L)
    return @. (A * spect.data - L*sqrt.(spect.data))/(4.0*pi)
end

function weyl_law(spect::DataSample, A, L, angles)
    return weyl_law(spect, A, L) .+ sum([(pi - phi/pi)/(24.0*phi) for phi in angles])
end

function weyl_law_k(spect::DataSample, A, L)
    return @. (A * spect.data^2 - L*spect.data)/(4.0*pi)
end

function weyl_law_k(spect::DataSample, A, L, angles)
    return weyl_law_k(spect, A, L) .+ sum([(pi - phi/pi)/(24.0*phi) for phi in angles])
end