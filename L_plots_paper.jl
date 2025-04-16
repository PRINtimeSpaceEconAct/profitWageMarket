set_theme!(Theme(
    fontsize = 20, # Set default font size for the whole figure
    palette = (
        color = [:black])  # all elements will be black
))

# Set it globally
set_theme!(bw_theme)
function plot1(sol,p,T1,T2)

    @unpack_parameters p
    profits1 = computeProfits(sol(T1),p)
    profits2 = computeProfits(sol(T2),p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    profitsEq = computeProfits(muEq,p)
    
    
    muMin = 0.9*min(minimum(sol(T1)),minimum(sol(T2)),minimum(muEq))
    muMax = 1.1*max(maximum(sol(T1)),maximum(sol(T2)),minimum(muEq))
    profitsMin = 0.9*min(minimum(profits1),minimum(profits2),minimum(profitsEq))
    profitsMax = 1.1*max(maximum(profits1),maximum(profits2),maximum(profitsEq))

    # Create a figure and axes
    fig = Figure(size = (2000, 1000))
    ax = Axis(fig[1, 1], xlabel = "Type of goods", ylabel = "Firms density μ", title = "",limits = (0, n, muMin, muMax))
    ax2 = Axis(fig[1, 2], xlabel = "Type of goods", ylabel = "Profits π", title = "",limits = (0, n, profitsMin, profitsMax))


    # Create lines for each plot
    mu1_line = lines!(ax, x, sol(T1), linewidth = 3, label = "Firms' density at t = " * string(round(T1,digits = 3)),linestyle = :solid)
    mu2_line = lines!(ax, x, sol(T2), linewidth = 3, label = "Firms' density at t = " * string(round(T2,digits = 3)),linestyle = :dash)
    profits1_line = lines!(ax2, x, profits1, linewidth = 3, label = "Profits per sector at t = " * string(round(T1,digits = 3)),linestyle = :solid)
    profits2_line = lines!(ax2, x, profits2, linewidth = 3, label = "Profits per sector at t = " * string(round(T2,digits = 3)),linestyle = :dash)
    firm_eq = lines!(ax, x, muEq, linewidth = 3, label = "Firms' density equilibirum",linestyle = :dot)
    profit_eq = lines!(ax2, x, profitsEq, linewidth = 3, label = "Firms' profits equilibirum",linestyle = :dot)

    # Add legend
    axislegend(ax)
    axislegend(ax2)

    # Display figure
    display(fig)

    Makie.save("fig1.png",fig)
    return nothing
end

function plot2(sol,p,T1,T2,T3)

    @unpack_parameters p
    profits0 = computeProfits(sol(0),p)
    profits1 = computeProfits(sol(T1),p)
    profits2 = computeProfits(sol(T2),p)
    profits3 = computeProfits(sol(T3),p)
    Δprofits1 = profits1 - profits0
    Δprofits2 = profits2 - profits0
    Δprofits3 = profits3 - profits0

    profits0Min = minimum(profits0)
    profits0Max = maximum(profits0)
    ΔprofitsMin = min(minimum(Δprofits1),minimum(Δprofits2),minimum(Δprofits3))
    ΔprofitsMax = max(maximum(Δprofits1),maximum(Δprofits2),maximum(Δprofits3))
    
    # Create a figure and axes
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], xlabel = "Profits π", ylabel = "Variation of profits Δπ",
         title = "",limits = (profits0Min, profits0Max, ΔprofitsMin, ΔprofitsMax))
    

    profits1_line = scatter!(ax, profits0, Δprofits1, label = "π(" * string(round(T1,digits=3)) * ") - π(0)", marker = :circle)
    profits2_line = scatter!(ax, profits0, Δprofits2, label = "π(" * string(round(T2,digits=3)) * ") - π(0)", marker = :star4)
    profits3_line = scatter!(ax, profits0, Δprofits3, label = "π(" * string(round(T3,digits=3)) * ") - π(0)", marker = :dtriangle)

    axislegend(ax)
    display(fig)

    Makie.save("fig2.png",fig)
end

function plot3()

    paramList = (parameters(), parameters(β = 0.8), parameters(σ = 0.3), parameters(η = 0.5))
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], xlabel = "Period t", ylabel = "log ||μ(t)-μEQ||₂", title = "")

    p = paramList[1]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    errNorm2 = [sqrt.(sum((sol(t) - muEq).^2*Δx)) for t in sol.t]
    lines!(ax, sol.t, log.(errNorm2), linewidth = 3, label = "Error norm baseline",linestyle = :solid)

    p = paramList[2]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    errNorm2 = [sqrt.(sum((sol(t) - muEq).^2*Δx)) for t in sol.t]
    lines!(ax, sol.t, log.(errNorm2), linewidth = 3, label = "Error norm β = 0.6",linestyle = :dash)

    p = paramList[3]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    errNorm2 = [sqrt.(sum((sol(t) - muEq).^2*Δx)) for t in sol.t]
    lines!(ax, sol.t, log.(errNorm2), linewidth = 3, label = "Error norm σ = 0.7",linestyle = :dot)

    p = paramList[4]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    errNorm2 = [sqrt.(sum((sol(t) - muEq).^2*Δx)) for t in sol.t]
    lines!(ax, sol.t, log.(errNorm2), linewidth = 3, label = "Error norm η = 0.3",linestyle = :dashdot)

    axislegend(ax)
    display(fig)

    Makie.save("fig3.png",fig)
end

function plot4()

    paramList = (parameters(), parameters(β = 0.8), parameters(σ = 0.3), parameters(η = -0.5))
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], xlabel = "Period t", ylabel = "log(XEq/X(t))", title = "")

    p = paramList[1]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    XEq = computeX(muEq,p)
    Xt = [computeX(sol(t),p) for t in sol.t]
    logErr = log.(XEq./Xt)
    lines!(ax, sol.t, logErr, linewidth = 3, label = "Error norm baseline",linestyle = :solid)

    p = paramList[2]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    XEq = computeX(muEq,p)
    Xt = [computeX(sol(t),p) for t in sol.t]
    logErr = log.(XEq./Xt)
    lines!(ax, sol.t, logErr, linewidth = 3, label = "Error norm β = 0.8",linestyle = :dash)

    p = paramList[3]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    XEq = computeX(muEq,p)
    Xt = [computeX(sol(t),p) for t in sol.t]
    logErr = log.(XEq./Xt)
    lines!(ax, sol.t, logErr, linewidth = 3, label = "Error norm σ = 0.3",linestyle = :dot)

    p = paramList[4]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    XEq = computeX(muEq,p)
    Xt = [computeX(sol(t),p) for t in sol.t]
    logErr = log.(XEq./Xt)
    lines!(ax, sol.t, logErr, linewidth = 3, label = "Error norm η = -0.5",linestyle = :dashdot)

    axislegend(ax)
    display(fig)

    Makie.save("fig4.png",fig)
end


function plot5(p,T1,T2,T3)

    sol,p = solveModelFixedCost(p)
    @unpack_parameters p

    profits0 = computeProfits(sol(0),p)
    profits1 = computeProfits(sol(T1),p)
    profits2 = computeProfits(sol(T2),p)
    profits3 = computeProfits(sol(T3),p)
    Δprofits1 = profits1 - profits0
    Δprofits2 = profits2 - profits0
    Δprofits3 = profits3 - profits0

    profits0Min = minimum(profits0)
    profits0Max = maximum(profits0)
    ΔprofitsMin = min(minimum(Δprofits1),minimum(Δprofits2),minimum(Δprofits3))
    ΔprofitsMax = max(maximum(Δprofits1),maximum(Δprofits2),maximum(Δprofits3))
    
    # Create a figure and axes
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], xlabel = "Profits π", ylabel = "Variation of profits Δπ",
         title = "",limits = (profits0Min, profits0Max, ΔprofitsMin, ΔprofitsMax))
    

    profits1_line = scatter!(ax, profits0, Δprofits1, label = "π(" * string(round(T1,digits=3)) * ") - π(0)", marker = :circle)
    profits2_line = scatter!(ax, profits0, Δprofits2, label = "π(" * string(round(T2,digits=3)) * ") - π(0)", marker = :star4)
    profits3_line = scatter!(ax, profits0, Δprofits3, label = "π(" * string(round(T3,digits=3)) * ") - π(0)", marker = :dtriangle)
     
    axislegend(ax)
    display(fig)

    Makie.save("fig5.png",fig)
end

function plot6(p,T1,T2,T3)
    

    sol,p = solveModelFixedCost(p)
    @unpack_parameters p

    profits1 = computeProfits(sol(T1),p)
    profits2 = computeProfits(sol(T2),p)
    profits3 = computeProfits(sol(T3),p)
    
    profitsMin = 0.9*min(minimum(profits1),minimum(profits2),minimum(profits3))
    profitsMax = 1.1*max(maximum(profits1),maximum(profits2),maximum(profits3))

    # Create a figure and axes
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], xlabel = "Type of goods", ylabel = "Profits π", title = "",limits = (0, n, profitsMin, profitsMax))

    # Create lines for each plot
    profits1_line = lines!(ax, x, profits1, linewidth = 3, label = "Profits per sector at t = " * string(round(T1,digits = 3)),linestyle = :solid)
    profits2_line = lines!(ax, x, profits2, linewidth = 3, label = "Profits per sector at t = " * string(round(T2,digits = 3)),linestyle = :dash)
    profits3_line = lines!(ax, x, profits3, linewidth = 3, label = "Profits per sector at t = " * string(round(T3,digits = 3)),linestyle = :dot)

    # Add legend
    axislegend(ax)

    # Display figure
    display(fig)

    Makie.save("fig6.png",fig)
    return nothing


end


plot1(sol,p,0,0.2)
plot2(sol,p,0.01,0.5,2)
plot3()
plot4()
plot5(p,0.01,0.5,2)
plot6(p,0,0.2,2)
