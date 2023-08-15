export GSE

"""
GSE <: Model 

`GSE` is a concrete type used to represent the Gaussian Symplectic Ensemble model of random matrix theory. 

## Description
Based on the Bohigas-Giannoni-Schmit conjecture this model describes the spectral statistics (in the semiclassical limit) 
of chaotic systems, invariant under time-reversal symmetry and including spin 1/2 interactions.

## Attributes
This model has no attributes.

## API
The following spectral statistcs can be evaluated for this model:
- [`level_spacing_pdf`](@ref)
- [`level_spacing_cdf`](@ref)
- [`level_spacing_u`](@ref)
- [`number_variance`](@ref)
- [`rigidity`](@ref)
- [`spectral_form_factor`](@ref)
"""
struct GSE <: Model end


#first argument of all functions is the variable


function level_spacing_pdf(model::GSE, s; n=1)
    #Wigner surmises for higher orders Wen-Jia Rao
    beta = 4.0 #level repulsion
    a = 0.5*n*(n+1.0)*beta + n -1
    A =  (gamma(0.5*a + 1.0)/gamma(0.5*a + 0.5) )^2.0 
    C = (2.0*(gamma(0.5*a + 1.0))^(a + 1.0)) / ((gamma(0.5*a + 0.5))^(a + 2.0))
    return @. C * s^a * exp(-A * s^2.0)
end

function level_spacing_cdf(model::GSE, s; n=1)
    beta = 4.0 #level repulsion
    a = 0.5*n*(n+1.0)*beta + n -1
    A =  (gamma(0.5*a + 1.0)/gamma(0.5*a + 0.5) )^2.0 
    C = (2.0*(gamma(0.5*a + 1.0))^(a + 1.0)) / ((gamma(0.5*a + 0.5))^(a + 2.0))
    return @. 0.5*C*s^(1.0 + a)*(A*s^2.0)^(0.5*(-1.0 - a))*(gamma((1.0 + a)*0.5) - gamma((1.0 + a)*0.5, A*s^2.0))

end

function level_spacing_u(model::GSE, s)
    cdf = level_spacing_cdf(model, s)
    return @. (2.0 / pi) * acos(sqrt(1.0 - cdf))    
end

function gap_probability(model::GSE, s)
    return @. exp(-64.0*s^2/(9.0*pi))*(9.0*pi+16.0* s^2.0)/(9.0*pi)-s*erfc(8.0/(3.0*sqrt(pi))*s)
end

function number_variance(model::GSE, l) 
    arg = 2.0 * pi .* l
    si, ci = sinint.(arg), cosint.(arg) 
    nv_gue = @. 1.0/pi^2.0 * (log(arg) + eulergamma + 1.0 - cos(arg) - ci) + l * (1.0 - 2.0/pi * si)    
    #si2 = sinint.(2.0*pi*l)
    return @. 0.5 * nv_gue + si^2.0/(4.0*pi^2.0) 
end 

function rigidity(model::GSE, l)
    beta = 4.0
    nv = @. 2.0/(beta*pi^2.0)*(log(beta*pi*l) + eulergamma + 1.0)
    #nv = @. 1.0/(2.0*pi^2.0)*(log(beta*pi*l) + eulergamma + 1.0 + pi^2.0/8.0)#approximate formula O(1/(2*pi*l)) 
    return @. 0.5 * nv - 9.0/(4.0*beta*pi^2.0)
end


function spectral_form_factor(model::GSE, t)
    sff = @. 0.5 * abs.(t) - 0.25*abs.(t)* log(abs(abs(t) - 1.0))
    idx = t.> 2.0
    sff[idx] .= 1.0    
    return sff
end