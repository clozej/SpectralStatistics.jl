
using StatsBase

function number_variance(E, L)
    #algorithm by Tomaz Prosen
    Ave1 = 0.0
    Ave2 = 0.0
    j = 2 #index of first state in interval
    k = 2 #index of last state in interval
    x = E[1] #current energy
    largest_energy = E[end - Int(ceil(L)+10)] #to make sure interval does not go out of energy bounds
    while x < largest_energy
        #move index k
        while E[k] < x+L
            k = k+1    
        end

        d1 = E[j] - x
        d2 = E[k] - (x+L)
        cn = k - j #number of states in interval
        if d1 < d2
            x = E[j]
            s = d1
            j = j + 1
        else
            x = E[k] - L
            s = d2
            k = k + 1
        end
        Ave1 = Ave1 + s*cn
        Ave2 = Ave2 + s*cn^2
    end
    #println(E[k])
    #println(largest_energy)
    s = largest_energy - E[1]
    Ave1 = Ave1/s
    Ave2 = Ave2/s
    AveSig = (Ave2 - Ave1^2)
    return AveSig
end

function number_variance(spect::S, x::Vector{T}) where {S<:DataSample, T<:Number} 
    E = spect.data
    Ls = x
    sz = length(x)
    nvs = zeros(sz) 
    Threads.@threads for i in 1:sz
        nvs[i] = number_variance(E, Ls[i])        
    end
    return x, nvs
end


NV = number_variance