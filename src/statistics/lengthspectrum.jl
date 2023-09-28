

function length_spectrum(l::Real, rho_fluct::Vector{T}, ks::Vector{T}) where T<:Real
    l_spec = complex(0.0)
    for (rho, k)  in zip(rho_fluct, ks) 
        l_spec += rho * exp(im * k * l)
    end
    norm = length(rho_fluct)
    return abs(l_spec)/norm
end


function length_spectrum(ls::Vector{T}, rho_fluct::Vector{T}, ks::Vector{T}) where T<:Real
    ls = collect(LinRange(min,max, grid))
    l_spec = zeros(grid) 
    Threads.@threads for i in 1:grid
        l_spec[i] = length_spectrum(ls[i],rho_fluct,ks)        
    end
    return l_spec
end