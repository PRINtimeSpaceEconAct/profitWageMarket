# Use a global, monochrome theme with larger fonts and thicker lines
Makie.set_theme!(Theme(
    fontsize = 32,
    palette = (color = [:black],), # force all series to black by default
    Lines = (linewidth = 5,),
    Scatter = (markersize = 14, strokewidth = 0.8),  # was 1.5, thinner outlines
    Axis = (
        xlabelsize = 36, ylabelsize = 36, titlesize = 40,
        xticklabelsize = 32, yticklabelsize = 32
    ),
    Legend = (
        labelsize = 28, titlesize = 30,
        patchsize = (110, 18),       # longer line segment in legend
        patchlabelgap = 16           # a bit more spacing after the line
    )
))

# set_theme!(Theme(
#     fontsize = 20, # Set default font size for the whole figure
#     palette = (
#         color = [:black])  # all elements will be black
# ))

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

    # μ figure (fig1a)
    figA = Figure(size = (1000, 1000))
    axA = Axis(figA[1, 1], xlabel = "Space of goods", ylabel = "Firms density μ",
               title = "", limits = (0, n, 0.8, 1.2))

    lines!(axA, x, sol(T1), linewidth = 5, color = :black,
           label = "t = " * string(round(T1, digits = 3)), linestyle = :solid)
    lines!(axA, x, sol(T2), linewidth = 5, color = :black,
           label = "t = " * string(round(T2, digits = 3)), linestyle = :dash)
    lines!(axA, x, muEq, linewidth = 5, color = :black,
           label = "t = ∞", linestyle = :dot)

    axislegend(axA)
    display(figA)
    Makie.save("fig1a.png", figA)

    # π figure (fig1b)
    figB = Figure(size = (1000, 1000))
    axB = Axis(figB[1, 1], xlabel = "Space of goods", ylabel = "Profit rate π",
               title = "", limits = (0, n, 0.11, 0.18))

    lines!(axB, x, profits1, linewidth = 5, color = :black,
           label = "t = " * string(round(T1, digits = 3)), linestyle = :solid)
    lines!(axB, x, profits2, linewidth = 5, color = :black,
           label = "t = " * string(round(T2, digits = 3)), linestyle = :dash)
    lines!(axB, x, profitsEq, linewidth = 5, color = :black,
           label = "t = ∞", linestyle = :dot)

    axislegend(axB)
    display(figB)
    Makie.save("fig1b.png", figB)

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
    ax = Axis(fig[1, 1], xlabel = "Profit rate π", ylabel = "Variation of profit rate Δπ",
         title = "",limits = (profits0Min, profits0Max, ΔprofitsMin, ΔprofitsMax))
    

    # Larger, black markers; at t = 0.5 use a hollow circle
    profits1_line = scatter!(ax, profits0, Δprofits1,
        label = "t = " * string(round(T1,digits = 3)),
        marker = isapprox(T1, 0.5; atol = 1e-6) ? :circle : :circle,
        color = isapprox(T1, 0.5; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)

    profits2_line = scatter!(ax, profits0, Δprofits2,
        label = "t = " * string(round(T2,digits = 3)),
        marker = isapprox(T2, 0.5; atol = 1e-6) ? :circle : :star4,
        color = isapprox(T2, 0.5; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)

    profits3_line = scatter!(ax, profits0, Δprofits3,
        label = "t = " * string(round(T3,digits = 3)),
        marker = isapprox(T3, 0.5; atol = 1e-6) ? :circle : :dtriangle,
        color = isapprox(T3, 0.5; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)

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
    lines!(ax, sol.t, log.(errNorm2), linewidth = 5, color = :black, label = "Baseline", linestyle = :solid)

    p = paramList[2]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    errNorm2 = [sqrt.(sum((sol(t) - muEq).^2*Δx)) for t in sol.t]
    lines!(ax, sol.t, log.(errNorm2), linewidth = 5, color = :black, label = "β = 0.6", linestyle = :dash)

    p = paramList[3]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    errNorm2 = [sqrt.(sum((sol(t) - muEq).^2*Δx)) for t in sol.t]
    lines!(ax, sol.t, log.(errNorm2), linewidth = 5, color = :black, label = "σ = 0.7", linestyle = :dot)

    p = paramList[4]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    errNorm2 = [sqrt.(sum((sol(t) - muEq).^2*Δx)) for t in sol.t]
    lines!(ax, sol.t, log.(errNorm2), linewidth = 5, color = :black, label = "η = 0.3", linestyle = :dashdot)

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
    lines!(ax, sol.t, logErr, linewidth = 5, color = :black, label = "Baseline", linestyle = :solid)
    minY = minimum(logErr)
    tmax = maximum(sol.t)

    p = paramList[2]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    XEq = computeX(muEq,p)
    Xt = [computeX(sol(t),p) for t in sol.t]
    logErr = log.(XEq./Xt)
    lines!(ax, sol.t, logErr, linewidth = 5, color = :black, label = "β = 0.8", linestyle = :dash)
    minY = min(minY, minimum(logErr))

    p = paramList[3]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    XEq = computeX(muEq,p)
    Xt = [computeX(sol(t),p) for t in sol.t]
    logErr = log.(XEq./Xt)
    lines!(ax, sol.t, logErr, linewidth = 5, color = :black, label = "σ = 0.3", linestyle = :dot)
    minY = min(minY, minimum(logErr))

    p = paramList[4]
    @unpack_parameters p
    sol,p = solveModel(p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    XEq = computeX(muEq,p)
    Xt = [computeX(sol(t),p) for t in sol.t]
    logErr = log.(XEq./Xt)
    lines!(ax, sol.t, logErr, linewidth = 5, color = :black, label = "η = -0.5", linestyle = :dashdot)
    minY = min(minY, minimum(logErr))

    ax.limits = (0, tmax, -0.001, 0.001)

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
    ax = Axis(fig[1, 1], xlabel = "Profit rate π", ylabel = "Variation of profit rate Δπ",
         title = "",limits = (profits0Min, profits0Max, ΔprofitsMin, ΔprofitsMax))
    

    # Larger, black markers; at t = 0.5 use a hollow circle
    profits1_line = scatter!(ax, profits0, Δprofits1,
        label = "t = " * string(round(T1,digits = 3)),
        marker = isapprox(T1, 0.1; atol = 1e-6) ? :circle : :circle,
        color = isapprox(T1, 0.1; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)

    profits2_line = scatter!(ax, profits0, Δprofits2,
        label = "t = " * string(round(T2,digits = 3)),
        marker = isapprox(T2, 0.1; atol = 1e-6) ? :circle : :star4,
        color = isapprox(T2, 0.1; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)

    profits3_line = scatter!(ax, profits0, Δprofits3,
        label = "t = " * string(round(T3,digits = 3)),
        marker = isapprox(T3, 0.1; atol = 1e-6) ? :circle : :dtriangle,
        color = isapprox(T3, 0.1; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)
     
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
    ax = Axis(fig[1, 1], xlabel = "Space of goods", ylabel = "Profit rate π", title = "",limits = (0, n, profitsMin, profitsMax))

    # Create lines for each plot (thicker, black, keep distinct linestyles)
    profits1_line = lines!(ax, x, profits1, linewidth = 5, color = :black, label = "t = " * string(round(T1,digits = 3)), linestyle = :solid)
    profits2_line = lines!(ax, x, profits2, linewidth = 5, color = :black, label = "t = " * string(round(T2,digits = 3)), linestyle = :dash)
    profits3_line = lines!(ax, x, profits3, linewidth = 5, color = :black, label = "t = " * string(round(T3,digits = 3)), linestyle = :dot)

    # Add legend
    axislegend(ax)

    # Display figure
    display(fig)

    Makie.save("fig6.png",fig)
    return nothing


end


function plot7(sol,p,T1,T2)

    @unpack_parameters p
    LabourProductivity1 = computeProduction(sol(T1),p)./p.L
    LabourProductivity2 = computeProduction(sol(T2),p)./p.L
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    LabourProductivityEq = computeProduction(muEq,p)./p.L

    production1 = computeProduction(sol(T1),p)
    production2 = computeProduction(sol(T2),p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    productionEq = computeProduction(muEq,p)

    # Y/L figure (fig7a)
    figA = Figure(size = (1000, 1000))
    axA = Axis(figA[1, 1], xlabel = "Space of goods", ylabel = "Value of production Y",
               title = "")

    lines!(axA, x, production1, linewidth = 5, color = :black,
           label = "t = " * string(round(T1, digits = 3)), linestyle = :solid)
    lines!(axA, x, production2, linewidth = 5, color = :black,
           label = "t = " * string(round(T2, digits = 3)), linestyle = :dash)
    lines!(axA, x, productionEq, linewidth = 5, color = :black,
           label = "t = ∞", linestyle = :dot)

    axislegend(axA) 
    display(figA)
    Makie.save("fig7a.png", figA)

    # Y/L figure (fig7b)
    figB = Figure(size = (1000, 1000))
    axB = Axis(figB[1, 1], xlabel = "Space of goods", ylabel = "Labour productivity Y/L",
               title = "")

    lines!(axB, x, LabourProductivity1, linewidth = 5, color = :black,
           label = "t = " * string(round(T1, digits = 3)), linestyle = :solid)
    lines!(axB, x, LabourProductivity2, linewidth = 5, color = :black,
           label = "t = " * string(round(T2, digits = 3)), linestyle = :dash)
    lines!(axB, x, LabourProductivityEq, linewidth = 5, color = :black,
           label = "t = ∞", linestyle = :dot)

    axislegend(axB)
    display(figB)
    Makie.save("fig7b.png", figB)

    return nothing
end

function plot8(sol,p,T1,T2,T3)

    @unpack_parameters p
    LabourProductivity0 = computeProduction(sol(0),p)./p.L
    LabourProductivity1 = computeProduction(sol(T1),p)./p.L
    LabourProductivity2 = computeProduction(sol(T2),p)./p.L
    LabourProductivity3 = computeProduction(sol(T3),p)./p.L
    ΔLabourProductivity1 = log.(LabourProductivity1) - log.(LabourProductivity0)
    ΔLabourProductivity2 = log.(LabourProductivity2) - log.(LabourProductivity0)
    ΔLabourProductivity3 = log.(LabourProductivity3) - log.(LabourProductivity0)

    LabourProductivity0Min = minimum(LabourProductivity0)
    LabourProductivity0Max = maximum(LabourProductivity0)
    ΔLabourProductivityMin = min(minimum(ΔLabourProductivity1),minimum(ΔLabourProductivity2),minimum(ΔLabourProductivity3))
    ΔLabourProductivityMax = max(maximum(ΔLabourProductivity1),maximum(ΔLabourProductivity2),maximum(ΔLabourProductivity3))
    
    # Create a figure and axes
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], xlabel = "Labour productivity Y/L", ylabel = "Growth rate of labour productivity",
         title = "",limits = (LabourProductivity0Min, LabourProductivity0Max, ΔLabourProductivityMin, ΔLabourProductivityMax))
    

    # Larger, black markers; at t = 0.5 use a hollow circle
    profits1_line = scatter!(ax, LabourProductivity0, ΔLabourProductivity1,
        label = "t = " * string(round(T1,digits = 3)),
        marker = isapprox(T1, 0.1; atol = 1e-6) ? :circle : :circle,
        color = isapprox(T1, 0.1; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)

    profits2_line = scatter!(ax, LabourProductivity0, ΔLabourProductivity2,
        label = "t = " * string(round(T2,digits = 3)),
        marker = isapprox(T2, 0.1; atol = 1e-6) ? :circle : :star4,
        color = isapprox(T2, 0.1; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)

    profits3_line = scatter!(ax, LabourProductivity0, ΔLabourProductivity3,
        label = "t = " * string(round(T3,digits = 3)),
        marker = isapprox(T3, 0.1; atol = 1e-6) ? :circle : :dtriangle,
        color = isapprox(T3, 0.1; atol = 1e-6) ? :transparent : :black,
        strokecolor = :black,
        markersize = 14)

    axislegend(ax)
    display(fig)

    Makie.save("fig8.png",fig)
end


function plot9()

    p = parameters(T_end=5.0)
    @unpack_parameters p
    sol,p = solveModel(p)
    solMobile,p = solveModelMobile(p)


    Xt = [computeX(sol(t),p) for t in sol.t]
    XtMobile = [computeXMobile(solMobile(t),p) for t in solMobile.t]

    Diff = XtMobile .- Xt

    fig = Figure(size = (1000, 1000))
    # ax = Axis(fig[1, 1], xlabel = "Period t", ylabel = "Efficiency loss ΔX", title = "")
    ax = Axis(fig[1, 1], xlabel = "Period t", ylabel = "Efficiency loss ΔX", title = "")
    # lines!(ax, sol.t, Xt, linewidth = 5, color = :black, label = "Immobile labour", linestyle = :solid)
    # lines!(ax, solMobile.t, XtMobile, linewidth = 5, color = :black, label = "Mobile labour", linestyle = :dash)

    lines!(ax, sol.t, (Diff), linewidth = 5, color = :black, label = "Immobile labour", linestyle = :solid)

    # axislegend(ax)
    display(fig)

    Makie.save("fig9.png",fig)
end

# p = parameters()
# sol,p = solveModel(p)
# plot1(sol,p,0,0.2)
# plot2(sol,p,0.01,0.5,2)
# plot3()
# plot4()
# sol,p = solveModelFixedCost(p)
# plot5(p,0.01,0.1,2)
# plot6(p,0,0.2,2)
# sol,p = solveModel(p)
# plot7(sol,p,0,0.2)
# plot8(sol,p,0.01,0.5,2)
plot9()

