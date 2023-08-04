function split_spectrum(spect, cut_value::T) where T<:Real
    idx = searchsortedfirst(spect, cut_value)
    return [spect[1:idx-1], spect[idx:end]]
end

function split_spectrum(spect, cut_values::AbstractArray)
    s1, s2 = split_spectrum(spect, cut_values[1])
    s_array = [s1]
    for cut in cut_values[2:end]
        s1, s2 = split_spectrum(s2, cut)
        push!(s_array,s1)
    end
    push!(s_array,s2)
    return s_array
end

function match_split(a, split)
    indices = [collect(1:length(s)) for s in split]
    for i in 2:length(indices)
        idx = indices[i]
        j = indices[i-1][end] #last index 
        indices[i] = idx .+ j
        #println(indices[i][1],indices[i][end])
    end
    return [a[idx] for idx in indices] 
end
