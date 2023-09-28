export Gamma

"""
Gamma <: Model 

`Gamma` is a concrete type used to represent the Gamma model. 

## Description
This model represents the gamma distributions. 

## Attributes
* `gamma`: The level repulsion exponent.

## API
The following spectral statistcs can be evaluated for this model:
- [`level_spacing_pdf`](@ref)
- [`level_spacing_cdf`](@ref)
"""
struct Gamma <: Model 
    gamma::Float64
end
#Gamma distributions of Bogomolny and Giraud
Gamma(;gamma=1.0) = Gamma(gamma)
Gamma((d::Dict)) = Gamma(d[:gamma])

function level_spacing_pdf(model::Gamma, s)
    g = model.gamma
    A = g + 1.0
    B = g + 1.0
    C = B^A / gamma(A)
    return @. C * s^(A - 1.0) * exp(-B * s)
end

function level_spacing_cdf(model::Gamma, s)
    g = model.gamma
    A = g + 1.0
    B = g + 1.0
    C = 1.0 / gamma(A)
    return [C*gamma_inc(A, B * si)[1] for si in s]
end

function level_spacing_u(model::Gamma, s)
    cdf = level_spacing_cdf(model, s)
    return @. (2.0 / pi) * acos(sqrt(abs(1.0 - cdf)))
end
