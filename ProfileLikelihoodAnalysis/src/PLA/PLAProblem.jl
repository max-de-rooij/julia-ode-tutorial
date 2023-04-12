using Optimization, Distributions, OptimizationOptimJL, LineSearches, DataInterpolations
include("PLAResult.jl")

mutable struct PLAProblem{T}
    loss_function
    calibrated_parameters::T
end

function __create_pla_loss(loss_fnc, par_idx, hyperparameters)
    loss(p, x) = begin
        parameters = [p[1:par_idx-1]; x; p[par_idx:end]]
        loss_fnc(parameters, hyperparameters)
    end
end

function __pla_step(optf::OptimizationFunction, par_value, calibrated_parameters, lb, ub)
    opt_problem = Optimization.OptimizationProblem(optf, calibrated_parameters, par_value, lb = lb, ub=ub);
    Optimization.solve(opt_problem, BFGS(linesearch=BackTracking(order=3)))
end

function __compute_confidence_interval(result::PLAResult)
    interp = LinearInterpolation(result.p[end-1:end], result.obj[end-1:end])
    interp(result.threshold+result.initial_obj)
end


function __pla(optf, calibrated_parameter_set, pla_parameter, loss_fnc, min_par, max_par, threshold, lb, ub, min_par_step, max_par_step, min_obj_change, max_obj_change, sample_size, direction)
    
    result = PLAResult([pla_parameter], [loss_fnc(calibrated_parameter_set, pla_parameter)[1]], pla_parameter, loss_fnc(calibrated_parameter_set, pla_parameter)[1], threshold, (-Inf, Inf))
    count = 1
    step = min_par_step
    res = result.obj[1]
    reset_max_change = false

    continuation = (pla_parameter, count, result, res, direction) -> begin
        edges = direction == :left ? pla_parameter >= min_par : pla_parameter <= max_par
        edges && (count <= sample_size) && (previous_obj(result)-res) < threshold
    end

    while continuation(pla_parameter, count, result, res, direction)

        temp_par = direction == :left ? previous_par(result)-step : previous_par(result) + step
        sol = __pla_step(optf, temp_par, calibrated_parameter_set, lb, ub)

        obj_change = abs(sol.objective - previous_obj(result))

        if step == min_par_step && obj_change > max_obj_change
            println("maximum objective function change encountered with minimum step size.\nTemporarily raising max_obj_change from $(max_obj_change) to $(obj_change)")
            original_max_change = max_obj_change
            max_obj_change = obj_change
            reset_max_change = true
        end

        if obj_change > max_obj_change
            step /= 2
            step = max(step, min_par_step)
        else
            count += 1
            append!(result, temp_par, sol.objective)
            pla_parameter = temp_par
            calibrated_parameter_set = sol.u

            if obj_change < min_obj_change
                step *= 2
                step = min(step, max_par_step)
            end
        end

        if reset_max_change
            max_obj_change = original_max_change
            reset_max_change = false
        end
    end

    if previous_obj(result)-res > threshold
        result.confint = direction == :left ? (__compute_confidence_interval(result), Inf) : (-Inf, __compute_confidence_interval(result))
    end
    result
end


solve(prob::PLAProblem, par_idx; min_par = :auto, max_par = :auto, threshold = :auto, lb = :auto, ub = :auto, 
    min_par_step = :auto, max_par_step = :auto, min_obj_change = :auto, max_obj_change = :auto, sample_size = 100, hyperparameters=nothing) = begin
    

    # minimum and maximum parameter values 
    min_par = min_par == :auto ? 0. : min_par[par_idx]  # Choose the correct min_par if it is not auto-select 
    max_par = max_par == :auto ? 10. .* prob.calibrated_parameters[par_idx] : max_par[par_idx]

    # Re-adjust the indices of the lower and upper bounds (choosing a subset without lb/ub of par_idx)
    lb = (lb == :auto ? nothing : lb[1:end .!= par_idx])
    ub = (ub == :auto ? nothing : ub[1:end .!= par_idx])
    
    
    # threshold
    threshold = threshold == :auto ? cquantile(Chisq(length(prob.calibrated_parameters)), 0.95) : threshold
    
    # parameter steps
    min_par_step = min_par_step == :auto ? maximum([prob.calibrated_parameters[par_idx]-min_par, max_par-prob.calibrated_parameters[par_idx]])/sample_size : min_par_step
    max_par_step = max_par_step == :auto ? maximum([prob.calibrated_parameters[par_idx]-min_par, max_par-prob.calibrated_parameters[par_idx]])/(0.1*sample_size) : max_par_step

    # objective function changes
    min_obj_change = min_obj_change == :auto ? 0.001 : min_obj_change
    max_obj_change = max_obj_change == :auto ? 0.05 : max_obj_change

    # create loss and optimization functions
    loss_fnc = __create_pla_loss(prob.loss_function, par_idx, hyperparameters)
    adtype = Optimization.AutoZygote()
    optf   = Optimization.OptimizationFunction(loss_fnc, adtype)

    calibrated_parameter_set = [prob.calibrated_parameters[1:par_idx-1]; prob.calibrated_parameters[par_idx+1:end]]
    pla_parameter = prob.calibrated_parameters[par_idx]

    left_result = __pla(optf, calibrated_parameter_set, pla_parameter, loss_fnc, min_par, max_par, threshold, lb, ub, min_par_step, max_par_step, min_obj_change, max_obj_change, sample_size, :left)
    right_result =  __pla(optf, calibrated_parameter_set, pla_parameter, loss_fnc, min_par, max_par, threshold, lb, ub, min_par_step, max_par_step, min_obj_change, max_obj_change, sample_size, :right)

    PLAResult(
        [reverse(left_result.p); right_result.p[2:end]],
        [reverse(left_result.obj); right_result.obj[2:end]],
        pla_parameter,
        left_result.initial_obj,
        threshold,
        (left_result.confint[1], right_result.confint[2])
    )
end
