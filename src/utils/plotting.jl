using Makie

default_linewidth = 0.75


#simple line plots
function plot_stat!(ax, x, y; 
    logy = false, logx = false,
    lineargs = Dict(:linewidth=>default_linewidth)
    )
        
    if logx
        x = log10.(x)
    end    
    
    if logy
        y = log10.(y)
    end
    
    lines!(ax,x,y;lineargs... )
end


function plot_stat!(ax, spectrum::D, statistic; limits=(0.0, 5.0), grid=200,
    logy = false, logx = false,
    lineargs = Dict(:color=>:black, :linewidth=>default_linewidth),
    statargs = Dict()
    ) where {D<:DataSample}
    
    x = collect(range(limits[1], limits[2], grid))
    y = statistic(spectrum, x; statargs...)
    
    plot_stat!(ax,x,y;logx,logy,lineargs)
end

function plot_stat!(ax, model::M, statistic; limits=(0.0, 5.0), grid=200,
    logy = false, logx = false,
    lineargs = Dict(:linewidth=>default_linewidth),
    statargs = Dict()
    ) where {M <:Model}
    
    x = collect(range(limits[1], limits[2], grid))
    y = statistic(model, x; statargs...)
    
    plot_stat!(ax,x,y;logx,logy,lineargs)
end

#plot difference

function plot_diff!(ax, spectrum::D, model::M, statistic; 
    limits=(0.0, 5.0), grid=200,
    lineargs = Dict(:linewidth=>default_linewidth),
    statargs = Dict()) where {D<:DataSample, M <:Model}

    x = collect(range(limits[1], limits[2], grid))
    y1 = statistic(spectrum, x; statargs...)
    y2 = statistic(model, x; statargs...)
    y = y1 - y2
    lines!(ax,x,y;lineargs... )

end

function plot_diff!(ax, x, y1, model::M, statistic; 
    limits=(0.0, 5.0), grid=200,
    lineargs = Dict(:linewidth=>default_linewidth),
    statargs = Dict()) where {M <:Model}

    x = collect(range(limits[1], limits[2], grid))
    #y1 = statistic(spectrum, x; statargs...)
    y2 = statistic(model, x; statargs...)
    y = y1 - y2
    lines!(ax,x,y;lineargs... )

end