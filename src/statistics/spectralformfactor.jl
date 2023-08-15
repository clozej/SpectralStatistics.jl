export spectral_form_factor
#=
function spectral_form_factor(t::Float64, E::Vector{Float64})
    sff_imag = 0.0
    sff_real = 0.0
    for e in E
        sff_imag += sin(e * t * 2*pi)
        sff_real += cos(e * t * 2*pi)
    end

    return sff_real^2 + sff_imag^2
end
=#

#this version is slightly faster
function spectral_form_factor(E::Vector{T}, t::T) where T<:Real
    sff = complex(0.0)
    for e in E
        sff += exp(im * e * t * 2*pi)
    end
    return abs2(sff)
end

function spectral_form_factor(E::Vector{T}, t::T, f::Function; fkwargs::D) where {T<:Real, D<:Dict}
    sff = complex(0.0)
    for e in E
        sff += f(e;fkwargs...)*exp(im * e * t * 2*pi)
    end
    return abs2(sff)
end


function spectral_form_factor(spect::UnfoldedSpectrum, ts::Vector{T}) where T<:Real
    E =spect.data
    grid = length(ts)
    sff = zeros(grid) 
    Threads.@threads for i in 1:grid
        sff[i] = spectral_form_factor(E, ts[i])        
    end
    return ts, sff
end

SFF = spectral_form_factor