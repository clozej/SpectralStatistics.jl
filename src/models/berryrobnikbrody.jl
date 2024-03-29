
export BerryRobnikBrody

"""
BerryRobnikBrody <: Model 

`BerryRobnikBrody` is a concrete type used to represent the Berry-Robnik-Brody model. 

## Description
This model is the Berry-Robnik model (see [`BerryRobnik`](@ref)) where instead of the GOE the [`Brody`](@ref) model is used for the chaotic part.   
It is commonly used to describe systems with a divided phase space as well as some degree of localization. 

## Attributes
* `rho`: The Liouville measure of the combined regular component.
* `beta`: The level repulsion exponent.

## API
The following spectral statistcs can be evaluated for this model:
- [`level_spacing_pdf`](@ref)
- [`level_spacing_cdf`](@ref)
- [`level_spacing_u`](@ref)
"""
struct BerryRobnikBrody <: Model 
    rho::Float64
    beta::Float64
end


BerryRobnikBrody(;rho=0.0,beta=1.0) = BerryRobnikBrody(rho, beta)
BerryRobnikBrody(d::Dict) = BerryRobnikBrody(d[:rho], d[:beta])
#first argument of all functions is the model
#definitions by Benjamin Batistić 
# we use the Brody distributions for the chaotic contributions

function level_spacing_pdf(model::BerryRobnikBrody, s)
    # rho is Liouville measure of the regular component
    # beta is the level repulsion parameter
    rho = model.rho
    beta = model.beta
    brody = Brody(beta)
    arg = @. (1.0-rho)*s
    gap = gap_probability(brody,arg)
    pdf = level_spacing_pdf(brody,arg)
    cdf = level_spacing_cdf(brody,arg)
    a = @. (rho^2.0)*gap-2.0*rho*(1.0-rho)*(cdf - 1.0) + ((1 - rho)^2.0)*pdf 
    return @. a*exp(-rho*s)

end

function level_spacing_cdf(model::BerryRobnikBrody, s)
    rho = model.rho
    beta = model.beta
    brody = Brody(beta)
    arg = @. (1.0-rho)*s
    gap = gap_probability(brody,arg)
    pdf = level_spacing_pdf(brody,arg)
    cdf = level_spacing_cdf(brody,arg)
    a = @. (1.0-rho)*(cdf-1.0)-rho*gap
    return @. a*exp(-rho*s) + 1.0
end

function level_spacing_u(model::BerryRobnikBrody, s)
    rho = model.rho
    beta = model.beta
    cdf = level_spacing_cdf(model, s)
    return @. (2.0 / pi) * acos(sqrt(abs(1.0 - cdf)))        
end

function gap_probability(model::BerryRobnikBrody, s)
    rho = model.rho
    beta = model.beta
    brody = Brody(beta)
    gap_chaotic = gap_probability(brody, (1.0 - rho)*s)
    return @. gap_chaotic*exp(-s*rho)
end