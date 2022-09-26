using StatsBase

#Ideally one would do this programatically using metaprograming

function pdf_hist(var, bins)
    h = StatsBase.fit(Histogram, var, bins)
    h = StatsBase.normalize(h, mode=:pdf)
    return midpoints(bins), h.weights
end

function cdf(var, pts)
    f = StatsBase.ecdf(var)
    return pts, f(pts)
end

function u_cdf(var, pts)
    f = StatsBase.ecdf(var)
    U(w) = (2.0 / pi) * arccos.(sqrt.(1.0 - f(w)))
    return pts, U(pts)
end


#=
function staircase(x)
    step = [i for i in 1:length(x)]
    return x, step
end

function mode_fluctuations(x)
    y, step = staircase(x)
    return x, step .- 0.5 - y 
end

=#