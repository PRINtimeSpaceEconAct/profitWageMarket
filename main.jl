
include("L_LoadAll.jl")

function solveModel(p)
    T_span = (0.0,p.T_end)
    prob = ODEProblem(df!, p.μ₀, T_span, p)
    sol = solve(prob,Rosenbrock23(),saveat=p.t_save)
    
    return sol,p
end

function solveModelFixedCost(p)
    T_span = (0.0,p.T_end)
    prob = ODEProblem(dfFixedCost!, p.μ₀, T_span, p)
    sol = solve(prob,Rosenbrock23(),saveat=p.t_save)
    return sol,p
end

