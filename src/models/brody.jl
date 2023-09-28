
export Brody

"""
Brody <: Model 

`Brody` is a concrete type used to represent the Brody model. 

## Description
This model interpolates between the [`Poisson`](@ref) (``\\beta=0``) and Wigner-Dyson (``\\beta=1``) distributions. 
It is commonly used to describe systems with some degree of localization. 

## Attributes
* `beta`: The level repulsion exponent.

## API
The following spectral statistcs can be evaluated for this model:
- [`level_spacing_pdf`](@ref)
- [`level_spacing_cdf`](@ref)
- [`level_spacing_u`](@ref)
"""
struct Brody <: Model 
beta::Float64
end

Brody(;beta=1.0) = Brody(beta)
Brody(d::Dict) = Brody(d[:beta])
#first argument of all functions is the model 

function level_spacing_pdf(model::Brody, s)
    #Wigner surmises for higher orders Wen-Jia Rao
    beta = model.beta
    a = gamma((beta + 2.0) / (beta + 1.0))^(beta + 1.0)
    return @. (beta + 1.0) * a * s^beta * exp(-a * s^(beta + 1.0))
end

function level_spacing_cdf(model::Brody, s)
    beta = model.beta
    a = gamma((beta + 2.0) / (beta + 1.0))^(beta + 1.0)
    return @. 1.0-exp(-a*s^(beta+1.0))
end

function level_spacing_u(model::Brody, s)
    beta = model.beta
    cdf = level_spacing_cdf(model, s; beta)
    return @. (2.0 / pi) * acos(sqrt(abs(1.0 - cdf)))    
end

function gap_probability(model::Brody, s)
    beta = model.beta
    a = gamma((beta + 2.0)/(beta + 1.0))^(beta + 1.0)
    x = @. a*s^(beta + 1.0)
    gap = [ gamma_inc(1.0/(beta + 1.0), var)[1] for var in x]
    return  gap
end