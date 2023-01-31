import Base.append!
using Plots

mutable struct PLAResult{T}
    p::Vector{T}
    obj::Vector{T}
    initial_p::T
    initial_obj::T
    threshold::T
    confint::Tuple{T, T}
end

previous_par(x::PLAResult) = x.p[end]
previous_obj(x::PLAResult) = x.obj[end]

append!(x::PLAResult, par, obj) = begin
    Base.append!(x.p, par)
    Base.append!(x.obj, obj)
end

Plots.plot(result::PLAResult; xlabel="Parameter Value", ylabel="Objective Function Value", titlefontsize=10, titlefontfamily="Helvetica Bold", legend = true, title="", kwargs...) = begin
    pl = Plots.plot(result.p, result.obj, linewidth=2, labels="", xlabel=xlabel, ylabel=ylabel, palette=:seaborn_colorblind, title=title, titlefontsize=titlefontsize, titlefontfamily=titlefontfamily, kwargs...)
    Plots.scatter!(pl, [result.initial_p], [result.initial_obj], labels=legend ? "Optimized Parameter Value" : "", markershape=:star, markersize=8, palette=:seaborn_colorblind)
    Plots.hline!(pl, [result.initial_obj+result.threshold], labels=legend ? "Chi-square Threshold" : "", palette=:seaborn_colorblind)
    pl
end