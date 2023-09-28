export BerryRobnik

"""
    BerryRobnik <: Model 

`BerryRobnik` is a concrete type used to represent the Berry-Robnik model (with two components). 

## Description
This model describes the spectral statistics (in the semiclassical limit) 
of systems whose classical phase space features both regular and chaotic motion.
These systems are commonly refered to as systems with divided phase space or mixed-type systems.
The model rests on the argument that in the semiclassical limit the spectrum subspectra belonging to states that describe 
belonging to the states to the distinct types of motion (regular and chaotic) will decompose into separate components.
The regular part of the spectrum is modeled by [`Poisson`](@ref) statistics and the chaotic part by [`GOE`](@ref) statistics.

## Attributes
* `rho`: The Liouville measure of the combined regular component.

## API
The following spectral statistcs can be evaluated for this model:
- [`level_spacing_pdf`](@ref)
- [`level_spacing_cdf`](@ref)
- [`level_spacing_u`](@ref)
- [`number_variance`](@ref)
"""
struct BerryRobnik <: Model 
    # rho is Liouville measure of the regular component
    rho::Float64
end


BerryRobnik(;rho=0.0) = BerryRobnik(rho)
BerryRobnik(d::Dict) = BerryRobnik(d[:rho])
#BerryRobnik(d::Dict) = BerryRobnik(d["rho"])
#first argument of all functions is the model
#definitions by Benjamin Batistić 
# we use the Brody distributions at beta = 1.0 for the GOE contributions

function level_spacing_pdf(model::BerryRobnik, s)
    # rho is Liouville measure of the regular component
    # beta is the level repulsion parameter
    rho = model.rho
    beta = 1.0
    brody = Brody(beta)
    arg = @. (1.0-rho)*s
    gap = gap_probability(brody,arg)
    pdf = level_spacing_pdf(brody,arg)
    cdf = level_spacing_cdf(brody,arg)
    a = @. (rho^2.0)*gap-2.0*rho*(1.0-rho)*(cdf - 1.0) + ((1 - rho)^2.0)*pdf 
    return @. a*exp(-rho*s)

end

function level_spacing_cdf(model::BerryRobnik, s)
    rho = model.rho
    beta = 1.0
    brody = Brody(beta)
    arg = @. (1.0-rho)*s
    gap = gap_probability(brody,arg)
    pdf = level_spacing_pdf(brody,arg)
    cdf = level_spacing_cdf(brody,arg)
    a = @. (1.0 - rho)*(cdf - 1.0) - rho * gap
    return @. a*exp(-rho*s) + 1.0
end

function level_spacing_u(model::BerryRobnik, s)
    cdf = level_spacing_cdf(model, s)
    return @. (2.0 / pi) * acos(sqrt(abs(1.0 - cdf)))
end

function gap_probability(model::BerryRobnik, s)
    rho = model.rho
    beta = 1.0
    brody = Brody(beta)
    gap_chaotic = gap_probability(brody, (1.0 - rho)*s)
    return @. gap_chaotic*exp(-s*rho)
end

function number_variance(model::BerryRobnik, l)
    rho = model.rho
    goe = GOE()
    rho_cha = 1.0 - rho
    return @. l*rho + number_variance(goe, l*rho_cha)
end