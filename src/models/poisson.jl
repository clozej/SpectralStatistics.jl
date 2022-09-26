using SpecialFunctions 
using Base.MathConstants


struct Poisson <: Model end


#first argument of all functions is the variable


function level_spacing_pdf(model::Poisson, s ; n=1)
    return @. n^n/factorial(n-1) * s^(n-1) * exp(-n*s)
end

function level_spacing_cdf(model::Poisson, s; n=1)
    return @. 1.0 - gamma(n, n*s)/gamma(n)
end

function level_spacing_u(model::Poisson, s)
    cdf = level_spacing_cdf(model, s)
    return @. (2.0 / pi) * acos(sqrt(1.0 - cdf))
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