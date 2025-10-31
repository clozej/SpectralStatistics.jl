export Gamma, GammaLinear, fit_gamma_linear

function gamma_pdf(x, a, b)
    c = b^a / gamma(a)
    return @. c * x^(a - 1.0) * exp(-b * x)
end

function gamma_cdf(x, a, b)
    c = 1.0 / gamma(a)
    return [c*gamma_inc(a, b * xi)[1] for xi in x]
end

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

function level_spacing_pdf(model::Gamma, s; n::Int = 0)
    g = model.gamma
    a = g + 1.0
    b = (g + 1.0)/(n+1)
    return gamma_pdf(s, a, b)
end

function level_spacing_cdf(model::Gamma, s; n::Int = 0)
    g = model.gamma
    a = g + 1.0
    b = (g + 1.0)/(n+1)
    return gamma_cdf(s, a, b)
end

function level_spacing_u(model::Gamma, s)
    cdf = level_spacing_cdf(model, s)
    return @. (2.0 / pi) * acos(sqrt(abs(1.0 - cdf)))
end


"""
GammaLinear <: Model 

`GammaLinear` is a concrete type used to represent the GammaLinear model. 

## Description
This model represents the gamma distributions with linear scaling of shape perameters. 

## Attributes
* `gamma`: The level repulsion exponent of the nearest neighbour level spacings.
* `p`: Slope (proportionality parameter) of the higher order level spacings.

## API
The following spectral statistcs can be evaluated for this model:
- [`level_spacing_pdf`](@ref)
- [`level_spacing_cdf`](@ref)
- [`spectral_form factor`](@ref)
"""
struct GammaLinear <: Model 
    gamma::Float64
    p::Float64
end

GammaLinear(;gamma=1.0,p=2.0) = GammaLinear(gamma,p)
GammaLinear((d::Dict)) = GammaLinear(d[:gamma],d[:p])


function level_spacing_pdf(model::GammaLinear, s; n::Int = 0)
    g = model.gamma + model.p*n #linear dependence on
    a = g + 1.0
    b = (g + 1.0)/(n+1)
    return gamma_pdf(s, a, b)
end

function level_spacing_cdf(model::GammaLinear, s; n::Int = 0)
    g = model.gamma + model.p*n 
    a = g + 1.0
    b = (g + 1.0)/(n+1)
    return gamma_cdf(s, a, b)
end

function level_spacing_u(model::GammaLinear, s)
    cdf = level_spacing_cdf(model, s)
    return @. (2.0 / pi) * acos(sqrt(abs(1.0 - cdf)))
end

#formula by Bogomolny and Giraud
function SFF_formula_gamma(t, p, k)    
    function g(t) #laplace transform of gamma distributions
        return @. ((1.0 + t/p)^(p-k-1.0))/((1.0+t/p)^p -1.0)*exp(-t*(p-k-1.0)/(p+t)) 
    end
    return @. 1.0 + 2*real(g(2*pi*im*t)) 
end

function spectral_form_factor(model::GammaLinear, t)
    return SFF_formula_gamma(t, model.p, model.gamma) 
end

function fit_gamma_linear(spectrum; limits = (0.0, 12.0), nmax=5, grid = 200)
    ga = Gamma()
    par = Vector{Float64}(undef, 0)
    for n in 0:nmax
        x = collect(range(limits[1], limits[2], grid))
        y = level_spacing_pdf(spectrum, x; n=n)
        ga = fit_model(ga,level_spacing_pdf,  x, y; statargs=Dict(:n=>n) )
        println(ga.gamma)
        push!(par, ga.gamma)
    end
    ns = [float(i) for i in 0:nmax]
    p0 = [1.0, 0.0]
    model(x,p) = @. p[1]*x + p[2]
    fit = LsqFit.curve_fit(model, ns, par, p0)
    p, gamma = fit.param
    return GammaLinear(gamma,p)
end