using LsqFit
using InteractiveUtils

function fit_model(spectrum, model, statistic::T, var) where T<: Function
    x, y = statistic(spectrum, var)

    type = typeof(model)
    me = InteractiveUtils.methodswith(type, statistic)[1] #find method corresponding to model
    par_symbols = Base.kwarg_decl(me) #keyword argument symbols
    
    function f(s, fit_params) #fitting function
        dict = Dict([(key,val) for (key,val) in zip(par_symbols,fit_params)])
        return statistic(model, s; dict...)
    end

    p0 = [0.5 for i in 1:length(par_symbols)] #initial fitting parameters
    fit = LsqFit.curve_fit(f, x, y, p0)
    
    return Dict([(key,val) for (key,val) in zip(par_symbols,fit.param)])
end


function fit_model(model, statistic::T, x, y) where T<: Function
    type = typeof(model)
    me = InteractiveUtils.methodswith(type, statistic)[1] #find method corresponding to model
    par_symbols = Base.kwarg_decl(me) #keyword argument symbols
    
    function f(s, fit_params) #fitting function
        dict = Dict([(key,val) for (key,val) in zip(par_symbols,fit_params)])
        return statistic(model, s; dict...)
    end

    p0 = [0.5 for i in 1:length(par_symbols)] #initial fitting parameters
    fit = LsqFit.curve_fit(f, x, y, p0)
    
    return Dict([(key,val) for (key,val) in zip(par_symbols,fit.param)])
end