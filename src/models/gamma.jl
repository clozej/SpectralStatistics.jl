using SpecialFunctions 

struct Gamma <: Model end
#Gamma distributions of Bogomolny and Giraud

function level_spacing_pdf(model::Gamma, s; g)
    A = (g+1.0)
    C = A^(g + 1.0) / gamma(g + 1.0)
    return @. C * s^g * exp(-A * s)
end

function level_spacing_cdf(model::Gamma, s; g)
    A = (g+1.0)
    return @. gamma_inc(g, A*s)
end