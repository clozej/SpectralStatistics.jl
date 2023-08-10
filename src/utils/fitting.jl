using LsqFit
using InteractiveUtils

function fit_model(spectrum, model, statistic::F, var) where F<:Function
    x, y = statistic(spectrum, var)

    type = typeof(model)
    #me = InteractiveUtils.methodswith(type, statistic)[1] #find method corresponding to model
    par_symbols = fieldnames(type) #keyword argument symbols
    
    function f(s, fit_params) #fitting function
        dict = Dict([(key,val) for (key,val) in zip(par_symbols,fit_params)])
        #println(dict)
        m = type(dict)
        return statistic(m, s)
    end

    p0 = [0.5 for i in 1:length(par_symbols)] #initial fitting parameters
    fit = LsqFit.curve_fit(f, x, y, p0)
    params = Dict([(key,val) for (key,val) in zip(par_symbols,fit.param)])
    return type(params)
end


function fit_model(model, statistic::F, x, y) where F<:Function
    type = typeof(model)
    #me = InteractiveUtils.methodswith(type, statistic)[1] #find method corresponding to model
    par_symbols = fieldnames(type) #keyword argument symbols
    
    function f(s, fit_params) #fitting function
        dict = Dict([(key,val) for (key,val) in zip(par_symbols,fit_params)])
        #println(dict)
        m = type(dict)
        return statistic(m, s)
    end

    p0 = [0.5 for i in 1:length(par_symbols)] #initial fitting parameters
    fit = LsqFit.curve_fit(f, x, y, p0)
    params = Dict([(key,val) for (key,val) in zip(par_symbols,fit.param)])
    return type(params)
end