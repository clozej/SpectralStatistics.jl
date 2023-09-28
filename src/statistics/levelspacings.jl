export level_spacing, level_spacing_pdf, level_spacing_cdf, level_spacing_u, level_spacing_ratio

"""
    level_spacing(spect::UnfoldedSpectrum; n::Int = 1) → s::Vector

Return the level spacings of order `n`. 
The nearest neighbour level spacings are the default, given by n=1.

## Description
The level spacings of order n are given by the difference between the n-th consecutive levels of the spectrum. 
The nearest neighbour level spacing is considered most commonly. 

## Arguments
* `spect`: The unfolded energy spectrum, given as an instance of type [`UnfoldedSpectrum`](@ref).

## Keyword arguments
*  `n=1` : The order of the level spacings.

## Returns
*  `s` : Vector of level spacings.

"""
function level_spacing(spect::UnfoldedSpectrum; n::Int = 1)
    e = spect.data
    s = e[1 + n : end] .- e[1 : end-n]
    if n > 1
        s = s ./ float(n) 
    end
    return s
end

# make a macro that adds pdf to statitic name and computes it!
"""
    level_spacing_pdf(spect::UnfoldedSpectrum, bins::Vector; n::Int = 1) → s::Vector p::Vector

Return a histogram of the probability density function of the level spacings of order `n`.
The nearest neighbour level spacings are the default, given by n=1.

## Description
The nearest neighbour level spacing distributions are the most commonly studied spectral statistics.

## Arguments
* `spect`: The unfolded energy spectrum, given as an instance of type [`UnfoldedSpectrum`](@ref).

* `bins`: The boundaries of the bin positions.

## Keyword arguments
*  `n=1` : The order of the level spacings.

## Returns

*  `p` : Vector of the probability contained in each bin.
"""
function level_spacing_pdf(spect::UnfoldedSpectrum, pts::Vector{T}; n::Int = 1) where T<:Real
    s = level_spacing(spect; n=n)
    return pdf_hist(s, pts)
end

"""
    level_spacing_cdf(spect::UnfoldedSpectrum, pts::Vector; n::Int = 1) → s::Vector w::Vector

Return the cumulative density function of the level spacings of order `n` evaluated at positions `pts`.
The nearest neighbour level spacings are the default, given by n=1.

## Arguments
* `spect`: The unfolded energy spectrum, given as an instance of type [`UnfoldedSpectrum`](@ref).

* `pts`: The positions where cumulative density function should be evaluated.

## Keyword arguments
*  `n=1` : The order of the level spacings.

## Returns

*  `w` : Vector of the cumulative probabilities.
"""
function level_spacing_cdf(spect::UnfoldedSpectrum, pts::Vector{T}; n::Int = 1) where T<:Real
    s = level_spacing(spect; n=n)
    return cdf(s, pts)
end


"""
    level_spacing_u(spect::UnfoldedSpectrum, pts::Vector; n::Int = 1) → s::Vector u::Vector

Return the spectraly normalized cumulative density function of the nearest neighbour level spacings evaluated at positions `pts`.

## Description
The nearest neighbour level spacing distributions are the most commonly studied spectral statistics.
To normalize the relative fluctuations it is useful to perform the following nonlinear transformation

```math
U(s) := \\frac{2}{\\pi}\\arccos\\sqrt{1-W(s)},
```

where ``W(s)`` is the cumulative level spacing distribution.


## Arguments
* `spect`: The unfolded energy spectrum, given as an instance of type [`UnfoldedSpectrum`](@ref).

* `pts`: The positions where cumulative density function should be evaluated.

## Keyword arguments
*  `n=1` : The order of the level spacings.

## Returns

*  `u` : Vector of the cumulative probabilities.
"""
function level_spacing_u(spect::UnfoldedSpectrum, pts::Vector{T}; n::Int = 1) where T<:Real
    s = level_spacing(spect; n=n)
    return u_cdf(s, pts)
end


function level_spacing_ratio(spect::DataSample; shift::Int=1, n::Int = 1)
    s = level_spacing(spect, n=n)
    shifted = circshift(s,-shift)
    r = s[1:end-shift] ./ shifted[1:end-shift]
    return r
end
