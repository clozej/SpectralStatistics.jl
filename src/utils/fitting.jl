export fit_model

"""
    fit_model(spectrum, model, statistic::Function; limits=(0.0, 5.0), grid=200,
    statargs=Dict())) → model

Return the model with the parameters adjusted to optimally fit the statistic computed from the spectrum.
The statistic is computed on a linear grid and then fitted with the apropriate analytical expression from the model.

## Arguments
* `spectrum`: The energy spectrum, given as an instance of subtype [`DataSample`](@ref) compatible with the statistic.

* `model`: The model we wish to fit, given as an instance of subtype [`Model`](@ref).

* `statistic`: The function used to compute the spectral statistic we wish to fit.


## Keyword arguments
*  `limits=(0.0, 5.0)` : Limitng values of the argumet of the statistic we wish to compute.

*  `grid=200` : Evaluation grid of the statistic.

## Returns

*  `model` : A new instance of the model with adjusted parameters.
"""
function fit_model(spectrum, model, statistic::F; limits=(0.0, 5.0), grid=200,
    statargs=Dict()) where F<:Function
    x = collect(range(limits[1], limits[2], grid))
    y = statistic(spectrum, x; statargs...)

    type = typeof(model)
    #me = InteractiveUtils.methodswith(type, statistic)[1] #find method corresponding to model
    par_symbols = fieldnames(type) #keyword argument symbols
    #check if fitting is needed
    if isempty(par_symbols)
        return model
    end
    
    function f(s, fit_params) #fitting function
        dict = Dict([(key,val) for (key,val) in zip(par_symbols,fit_params)])
        #println(dict)
        m = type(dict)
        return statistic(m, s;  statargs...)
    end

    p0 = [0.5 for i in 1:length(par_symbols)] #initial fitting parameters
    fit = LsqFit.curve_fit(f, x, y, p0)
    params = Dict([(key,val) for (key,val) in zip(par_symbols,fit.param)])
    return type(params)
end

function fit_model(spectrum, model, statistic::F, var; statargs=Dict()) where F<:Function
    x = var
    y = statistic(spectrum, var; statargs...)

    type = typeof(model)
    #me = InteractiveUtils.methodswith(type, statistic)[1] #find method corresponding to model
    par_symbols = fieldnames(type) #keyword argument symbols
    #check if fitting is needed
    if isempty(par_symbols)
        return model
    end
    
    function f(s, fit_params) #fitting function
        dict = Dict([(key,val) for (key,val) in zip(par_symbols,fit_params)])
        #println(dict)
        m = type(dict)
        return statistic(m, s; statargs...)
    end

    p0 = [0.5 for i in 1:length(par_symbols)] #initial fitting parameters
    fit = LsqFit.curve_fit(f, x, y, p0)
    params = Dict([(key,val) for (key,val) in zip(par_symbols,fit.param)])
    return type(params)
end


function fit_model(model, statistic::F, x, y; statargs=Dict()) where F<:Function
    type = typeof(model)
    #me = InteractiveUtils.methodswith(type, statistic)[1] #find method corresponding to model
    par_symbols = fieldnames(type) #keyword argument symbols
    #check if fitting is needed
    if isempty(par_symbols)
        return model
    end
    
    function f(s, fit_params) #fitting function
        dict = Dict([(key,val) for (key,val) in zip(par_symbols,fit_params)])
        #println(dict)
        m = type(dict)
        return statistic(m, s; statargs...)
    end

    p0 = [0.5 for i in 1:length(par_symbols)] #initial fitting parameters
    fit = LsqFit.curve_fit(f, x, y, p0)
    params = Dict([(key,val) for (key,val) in zip(par_symbols,fit.param)])
    return type(params)
end