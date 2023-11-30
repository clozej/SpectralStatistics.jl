

export unfold_spectrum, weyl_law

"""
    unfold_spectrum(spect::DataSample, f::Function) → unfolded::UnfoldedSpectrum

Return the unfolded spectrum of `spect` by using the function `f` 
as the smooth part of the integrated density of states.

## Arguments

* `spect`: Spectrum to be unfolded.
*  `f(x)` : Function modeling the smooth part of the integrated density of states, where the argument `x` is the energy.

## Returns
*  `unfolded` : The unfolded spectrum as an instance of type [`UnfoldedSpectrum`](@ref).
"""
function unfold_spectrum(spect::DataSample, f::Function) #f is the unfolding function
    unfolded = @. f(spect.data)
    return UnfoldedSpectrum(unfolded)
end

"""
    unfold_spectrum(spect::DataSample, n::Int) → unfolded::UnfoldedSpectrum

Return the unfolded spectrum of `spect` by fitting a polynomial of degree `n` 
to the integrated density of states.

## Arguments
* `n` : Degree of polynomial moddeling the smooth part of the integrated density of states.
"""
function unfold_spectrum(spect::DataSample, n::Int ) #f is the unfolding function
    coeffs = fit_integrated_density(spect, n)
    unfolded = [evalpoly(x, coeffs) for x in spect.data]
    return UnfoldedSpectrum(unfolded)
end

"""
    unfold_spectrum(spect::DataSample, n::Int, cut_values) → unfolded::UnfoldedSpectrum

Return the unfolded spectrum of `spect` by piecewise fitting polynomials of degree `n` 
to the integrated density of states.

## Arguments
* `cut_values` : Relative positions of the cuts between the spectral pieces.  
"""
function unfold_spectrum(spect::DataSample, n::Int, cut_values) 
    x = spect.data
    l = x[end] - x[1]
    x = (x .- x[1]) ./ l
    y = [float(i) for i in 1:length(x)]
    y = y ./ length(y)
    split_x = split_spectrum(x, cut_values)
    split_y = match_split(y, split_x)
    fits = [Polynomials.fit(xi,yi,n+1).coeffs for (xi,yi) in zip(split_x,split_y)]  
    unfold_pieces = []
    for i in 1:length(split_x)
        x_i = split_x[i]
        y_i = split_y[i]
        coeffs = fits[i]
        unfolded = [evalpoly(xi, coeffs) for xi in x_i]
        append!(unfold_pieces,unfolded)
    end
    l = length(x)
    return UnfoldedSpectrum(l .* unfold_pieces)
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
    fit = Polynomials.fit(x,y,n+1)
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