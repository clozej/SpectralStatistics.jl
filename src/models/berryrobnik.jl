include("brody.jl")
include("goe.jl")

struct BerryRobnik <: Model end

#first argument of all functions is the model
#definitions by Benjamin Batistić 
# we use the Brody distributions at beta = 1.0 for the GOE contributions

function level_spacing_pdf(model::BerryRobnik, s; rho)
    # rho is Liouville measure of the regular component
    # beta is the level repulsion parameter
    beta = 1.0
    brody = Brody()
    arg = @. (1.0-rho)*s
    gap = gap_probability(brody,arg;beta)
    pdf = level_spacing_pdf(brody,arg;beta)
    cdf = level_spacing_cdf(brody,arg;beta)
    a = @. (rho^2.0)*gap-2.0*rho*(1.0-rho)*(cdf - 1.0) + ((1 - rho)^2.0)*pdf 
    return @. a*exp(-rho*s)

end

function level_spacing_cdf(model::BerryRobnik, s; rho)
    beta = 1.0
    brody = Brody()
    arg = @. (1.0-rho)*s
    gap = gap_probability(brody,arg;beta)
    pdf = level_spacing_pdf(brody,arg;beta)
    cdf = level_spacing_cdf(brody,arg;beta)
    a = @. (1.0 - rho)*(cdf - 1.0) - rho * gap
    return @. a*exp(-rho*s) + 1.0
end

function level_spacing_u(model::BerryRobnik, s; rho)
    cdf = level_spacing_cdf(model, s; rho)
    return @. (2.0 / pi) * acos(sqrt(1.0 - cdf))
end

function gap_probability(model::BerryRobnik, s; rho)
    beta = 1.0
    brody = Brody()
    gap_chaotic = gap_probability(brody, (1.0 - rho)*s; beta)
    return @. gap_chaotic*exp(-s*rho)
end

function number_variance(model::BerryRobnik, l; rho)
    goe = GOE()
    rho_cha = 1.0 - rho
    return @. l*rho + number_variance(goe, l*rho_cha)
end