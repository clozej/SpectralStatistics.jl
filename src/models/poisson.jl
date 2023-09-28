export Poisson

"""
Poisson <: Model 

`Poisson` is a concrete type used to represent the Poisson model. 

## Description
The Poisson model applies to sequences of independent random variables. 
Based on the Berry-Tabor conjecture the spectral statistics (in the semiclassical limit) 
of integrable systems are described by this model.

## Attributes
This model has no attributes.

## API
The following spectral statistcs can be evaluated for this model:
- [`level_spacing_pdf`](@ref)
- [`level_spacing_cdf`](@ref)
- [`level_spacing_u`](@ref)
- [`number_variance`](@ref)
- [`rigidity`](@ref)
- [`spectral_form_factor`](@ref)
"""
struct Poisson <: Model end


#first argument of all functions is the variable

"""
    level_spacing_pdf(model::Model, pts::Vector; n::Int = 1) → s::Vector p::Vector

Return the analytical expression for the level spacing probability density function, corresponding to the chosen model, evaluated at positions `pts`.

## Description
The nearest neighbour level spacing distributions are the most commonly studied spectral statistics.

## Arguments
* `model`: The model, given as an instance of a concrete subtype of [`Model`](@ref).

* `pts`: The positions where probabability density function should be evaluated.

## Keyword arguments
*  `n=1` : The order of the level spacings.

## Returns
*  `s` : Vector of the evaluation points.

*  `p` : Vector of the probabilites.
"""
function level_spacing_pdf(model::Poisson, s ; n=1)
    return @. n^n/factorial(n-1) * s^(n-1) * exp(-n*s)
end

"""
    level_spacing_cdf(model::Model, pts::Vector; n::Int = 1) → s::Vector w::Vector

Return the analytical expression for level spacing cumulative density function, corresponding to the chosen model, evaluated at positions `pts`.

## Arguments
* `model`: The model, given as an instance of a concrete subtype of [`Model`](@ref).

* `pts`: The positions where cumulative density function should be evaluated.

## Keyword arguments
*  `n=1` : The order of the level spacings.

## Returns
*  `s` : Vector of the evaluation points.

*  `w` : Vector of the cumulative probabilities.
"""
function level_spacing_cdf(model::Poisson, s; n=1)
    return @. 1.0 - gamma(n, n*s)/gamma(n)
end

"""
    level_spacing_u(spect::UnfoldedSpectrum, pts::Vector; n::Int = 1) → s::Vector u::Vector

Return the analytical expression for the spectraly normalized cumulative density function of the nearest neighbour level spacings evaluated at positions `pts`.


## Arguments
* `model`: The model, given as an instance of a concrete subtype of [`Model`](@ref).

* `pts`: The positions where cumulative density function should be evaluated.

## Keyword arguments
*  `n=1` : The order of the level spacings.

## Returns
*  `s` : Vector of the evaluation points.

*  `u` : Vector of the cumulative probabilities.
"""
function level_spacing_u(model::Poisson, s)
    cdf = level_spacing_cdf(model, s)
    return @. (2.0 / pi) * acos(sqrt(abs(1.0 - cdf)))
end

function gap_probability(model::Poisson, s)
    return @. exp(-s)
end

function number_variance(model::Poisson, l) 
  return l
end 

function rigidity(model::Poisson, l) 
    return l ./ 15.0
  end 

function spectral_form_factor(model::Poisson, t)
    return ones(length(t))
end