using SpecialFunctions 

struct Gamma <: Model 
    gamma::Float64
end
#Gamma distributions of Bogomolny and Giraud
Gamma(;gamma=1.0) = Gamma(gamma)
Gamma((d::Dict)) = Gamma(d[:gamma])

function level_spacing_pdf(model::Gamma, s)
    g = model.gamma
    A = (g+1.0)
    C = A^(g + 1.0) / gamma(g + 1.0)
    return @. C * s^g * exp(-A * s)
end

function level_spacing_cdf(model::Gamma, s)
    g = model.gamma
    A = (g+1.0)
    return @. gamma_inc(g, A*s)
end

