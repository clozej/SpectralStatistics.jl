using SpecialFunctions 

struct Brody <: Model end


#first argument of all functions is the model 

function level_spacing_pdf(model::Brody, s; beta)
    #Wigner surmises for higher orders Wen-Jia Rao
    a = gamma((beta + 2.0) / (beta + 1.0))^(beta + 1.0)
    return @. (beta + 1.0) * a * s^beta * exp(-a * s^(beta + 1.0))
end

function level_spacing_cdf(model::Brody, s; beta)
    a = gamma((beta + 2.0) / (beta + 1.0))^(beta + 1.0)
    return @. 1.0-exp(-a*s^(beta+1.0))
end

function level_spacing_u(model::Brody, s; beta)
    cdf = level_spacing_cdf(model, s; beta)
    return @. (2.0 / pi) * acos(sqrt(1.0 - cdf))    
end

function gap_probability(model::Brody, s; beta)
    a = gamma((beta + 2.0)/(beta + 1.0))^(beta + 1.0)
    x = @. a*s^(beta + 1.0)
    gap = [ gamma_inc(1.0/(beta + 1.0), var)[1] for var in x]
    return  gap
end