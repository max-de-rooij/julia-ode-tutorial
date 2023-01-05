using DifferentialEquations, SciMLSensitivity, OptimizationOptimisers, Optimization, OptimizationOptimJL, LinearAlgebra, Zygote
using ForwardDiff
using DiffEqBase: UJacobianWrapper
using Plots
function rober(du,u,p,t)
    y₁,y₂,y₃ = u
    k₁,k₂,k₃ = p
    du[1] = -k₁*y₁+k₃*y₂*y₃
    du[2] =  k₁*y₁-k₂*y₂^2-k₃*y₂*y₃
    du[3] =  k₂*y₂^2
    nothing
end

p = [0.04,3e7,1e4]
u0 = [1.0,0.0,0.0]
prob = ODEProblem(rober,u0,(0.0,1e5),p)
sol = solve(prob,Rosenbrock23())
ts = sol.t
Js = map(u->I + 0.1*ForwardDiff.jacobian(UJacobianWrapper(rober, 0.0, p), u), sol.u)

function predict_adjoint(p)
    p = exp.(p)
    _prob = remake(prob,p=p)
    Array(solve(_prob,Rosenbrock23(autodiff=false),saveat=ts,sensealg=QuadratureAdjoint(autojacvec=ReverseDiffVJP(true))))
end

function loss_adjoint(p)
    prediction = predict_adjoint(p)
    prediction = [prediction[:, i] for i in axes(prediction, 2)]
    diff = map((J,u,data) -> J * (abs2.(u .- data)) , Js, prediction, sol.u)
    loss = sum(abs, sum(diff)) |> sqrt
    loss, prediction
end

callback = function (p,l,pred) #callback function to observe training
    println("Loss: $l")
    println("Parameters: $(exp.(p))")
    # using `remake` to re-create our `prob` with current parameters `p`
    plot(solve(remake(prob, p=exp.(p)), Rosenbrock23())) |> display
    return false # Tell it to not halt the optimization. If return true, then optimization stops
end

initp = ones(3)
# Display the ODE with the initial parameter values.
callback(initp,loss_adjoint(initp)...)

adtype = Optimization.AutoZygote()
optf = Optimization.OptimizationFunction((x,p)->loss_adjoint(x), adtype)
optprob = Optimization.OptimizationProblem(optf, initp)

# res = Optimization.solve(optprob, Optimisers.Adam(0.01), callback = callback, maxiters = 300)

# optprob2 = Optimization.OptimizationProblem(optf, res.u)

res2 = Optimization.solve(optprob, BFGS(), callback = callback, maxiters = 300, allow_f_increases=true)
println("Ground truth: $(p)\nFinal parameters: $(round.(exp.(res2.u), sigdigits=5))\nError: $(round(norm(exp.(res2.u) - p) ./ norm(p) .* 100, sigdigits=3))%")
