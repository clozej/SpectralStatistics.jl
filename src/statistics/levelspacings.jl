

function level_spacing(spect::DataSample; n::Int = 1)
    e = spect.data
    s = e[1 + n : end] .- e[1 : end-n]
    if n > 1
        s = s ./ float(n) 
    end
    return s
end

# make a macro that adds pdf to statitic name and computes it!
function level_spacing_pdf(spect::DataSample, bins::Vector{T}; n::Int = 1) where T
    s = level_spacing(spect; n=n)
    return pdf_hist(s, bins)
end

function level_spacing_cdf(spect::DataSample, pts::Vector{T}; n::Int = 1) where T
    s = level_spacing(spect; n=n)
    return cdf(s, pts)
end

function level_spacing_u(spect::DataSample, pts::Vector{T}; n::Int = 1) where T
    s = level_spacing(spect; n=n)
    return u_cdf(s, pts)
end


function level_spacing_ratio(spect::DataSample; shift::Int=1, n::Int = 1)
    s = level_spacing(spect, n=n)
    shifted = circshift(s,-shift)
    r = s[1:end-shift] ./ shifted[1:end-shift]
    return r
end


