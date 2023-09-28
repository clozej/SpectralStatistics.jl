using LsqFit
using InteractiveUtils


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