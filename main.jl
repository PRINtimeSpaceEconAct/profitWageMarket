include("L_LoadAll.jl")
using Pkg
Pkg.activate(@__DIR__)

"""
solveModel(p) -> (sol,p)
Baseline model solution.
"""
function solveModel(p)
    T_span = (0.0,p.T_end)
    prob = ODEProblem(df!, p.μ₀, T_span, p)
    sol = solve(prob,Rosenbrock23(),saveat=p.t_save)
    
    return sol,p
end

"""
solveModelMobile(p) -> (sol,p)
Mobile variant (uses dfMobile! and mobile profit / price index).
"""
function solveModelMobile(p)
    T_span = (0.0,p.T_end)
    prob = ODEProblem(dfMobile!, p.μ₀, T_span, p)
    sol = solve(prob,Rosenbrock23(),saveat=p.t_save)
    
    return sol,p
end


"""
solveModelFixedCost(p) -> (sol,p)
Fixed cost variant (uses dfFixedCost!).
"""
function solveModelFixedCost(p)
    T_span = (0.0,p.T_end)
    prob = ODEProblem(dfFixedCost!, p.μ₀, T_span, p)
    sol = solve(prob,Rosenbrock23(),saveat=p.t_save)
    return sol,p
end

"""
runAll(p)
Convenience: solves baseline, mobile, fixed–cost models and produces paper figures.
"""
function runAll(p)
    sol,p = solveModel(p)
    plot1(sol,p,0,0.2)
    plot2(sol,p,0.01,0.5,2)
    plot3()
    plot4()

    sol,p = solveModelFixedCost(p)
    plot5(p,0.01,0.1,2)
    plot6(p,0,0.2,2)

    sol,p = solveModel(p)
    plot7(sol,p,0,0.2)
    plot8(sol,p,0.01,0.5,2)
    plot9()
    
end