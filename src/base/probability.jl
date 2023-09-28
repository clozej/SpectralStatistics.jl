

#Ideally one would do this programatically using metaprograming

function pdf_hist(var, pts)
    bins = midpoints(pts)
    prepend!(bins, 0.0)
    last_dist = bins[end] - bins[end-1]
    last_pt = bins[end] + last_dist/2.0 
    append!(bins,last_pt)
    h = StatsBase.fit(Histogram, var, bins)
    h = StatsBase.normalize(h, mode=:pdf)
    return h.weights
end

function cdf(var, pts)
    f = StatsBase.ecdf(var)
    return f(pts)
end

function u_cdf(var, pts)
    f = StatsBase.ecdf(var)
    U(w) = @. (2.0 / pi) * acos(sqrt(abs(1.0 - f(w))))
    return U(pts)
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