include("brody.jl")
using SpecialFunctions 

struct BerryRobnikBrody <: Model end

#first argument of all functions is the model
#definitions by Benjamin Batistić 
# we use the Brody distributions for the chaotic contributions

function level_spacing_pdf(model::BerryRobnikBrody, s; rho, beta)
    # rho is Liouville measure of the regular component
    # beta is the level repulsion parameter
    brody = Brody()
    arg = @. (1.0-rho)*s
    gap = gap_probability(brody,arg;beta)
    pdf = level_spacing_pdf(brody,arg;beta)
    cdf = level_spacing_cdf(brody,arg;beta)
    a = @. (rho^2.0)*gap-2.0*rho*(1.0-rho)*(cdf - 1.0) + ((1 - rho)^2.0)*pdf 
    return @. a*exp(-rho*s)

end

function level_spacing_cdf(model::BerryRobnikBrody, s; rho, beta)
    brody = Brody()
    arg = @. (1.0-rho)*s
    gap = gap_probability(brody,arg;beta)
    pdf = level_spacing_pdf(brody,arg;beta)
    cdf = level_spacing_cdf(brody,arg;beta)
    a = @. (1.0-rho)*(cdf-1.0)-rho*gap
    return @. a*exp(-rho*s) + 1.0
end

function level_spacing_u(model::BerryRobnikBrody, s; rho, beta)
    cdf = level_spacing_cdf(model, s; rho, beta)
    return @. (2.0 / pi) * acos(sqrt(1.0 - cdf))        
end

function gap_probability(model::BerryRobnikBrody, s; rho, beta)
    brody = Brody()
    gap_chaotic = gap_probability(brody, (1.0 - rho)*s; beta)
    return @. gap_chaotic*exp(-s*rho)
end