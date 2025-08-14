"""
plotSol(sol, p; NSteps=1000, video=false, fName="animation.mp4")

- Plots the time evolution of firms' density μ(x, t) and profits π(x, t), alongside their steady states.
- If video=true, records an animation to fName.
"""
function plotSol(sol,p; NSteps = 1000, video = false, fName = "animation.mp4")

    @unpack_parameters p
    lmax = 1.1*max(maximum(sol))
    profits = computeProfits(sol(0),p)
    muEq = (A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) / sum((A .* L.^β).^((σ-1)/( σ*(β-η)*(σ-1)/σ + 1 ) ) * Δx)
    profitsEq = computeProfits(muEq,p)

    # Create a figure and axes
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], xlabel = "type of goods", ylabel = "", title = "Dynamic Plot",limits = (0, n, -0.5, lmax))

    # Create lines for each plot
    firm_density_line = lines!(ax, x, sol(0), linewidth = 3, label = "Firms' density")
    profits_line = lines!(ax, x, profits, linewidth = 3, label = "Profits per sector")
    derivative_profits_line = lines!(ax, x, ∂x(profits,p), linewidth = 3, label = "Derivative of profits per sector")
    firm_eq = lines!(ax, x, muEq, linewidth = 3, label = "Firms' density equilibirum")
    profit_eq = lines!(ax, x, profitsEq, linewidth = 3, label = "Firms' profits equilibirum")

    # Add legend
    axislegend(ax)

    # Display figure
    display(fig)


    # Loop through time steps
    for t in LinRange(0, p.T_end, NSteps)
        # Compute profits
        profits = computeProfits(sol(t),p)
        profits[isinf.(profits)] .= 0

        # Update the data in the lines
        firm_density_line[2] = sol(t)
        profits_line[2] = profits
        derivative_profits_line[2] = ∂x(profits,p)

        # Update title
        ax.title = "t = " * string(round(t, digits = 3))
        
        display(fig)
        sleep(0.01)
    end

    if video
        # Define the recording
        record(fig, fName, LinRange(0, p.T_end, NSteps)) do t
            # Compute profits
            profits = computeProfits(sol(t),p)
            profits[isinf.(profits)] .= 0
            # Update the data in the lines
            firm_density_line[2] = sol(t)
            profits_line[2] = profits
            derivative_profits_line[2] = ∂x(profits,p)

            # Update title
            ax.title = "t = " * string(round(t, digits = 3))
            
            display(fig)
            sleep(0.01)
        end
    end

    return nothing
end


"""
crossSectionalsPlot(sol, p; NSteps=1000, video=false, fName="crossPlot.mp4")

- Shows μ(x, t) and π(x, t) over x and their cross-sectional plots (sorted by value).
- If video=true, records an animation to fName.
"""
function crossSectionalsPlot(sol,p; NSteps = 1000, video = false, fName = "crossPlot.mp4")

    @unpack_parameters p
    # f(x) = -tan(5.5*(x-1/4))+30
    
    profits = computeProfits(sol(0),p)
    profits[isinf.(profits)] .= 0    
    crossμX0,crossμY = toCrossectional(sol(0),p)
    crossProfitsX0,crossProfitsY = toCrossectional(profits,p)

    distributionAreaProfits =   cumsum( crossProfitsY[2:end] .* diff(crossProfitsX0) )
    iLower = findmin(distributionAreaProfits .< 0.1)[2]
    iUpper = findmax(distributionAreaProfits .> 0.9)[2]
    crossProfitsX0[1:iLower] .= crossProfitsX0[iLower]
    crossProfitsX0[iUpper:end] .= crossProfitsX0[iUpper]
    crossProfitsY[1:iLower] .= crossProfitsY[iLower]
    crossProfitsY[iUpper:end] .= crossProfitsY[iUpper]

    # Create a figure and axes
    fig = Figure(size = (1000, 1000))

    ax = Axis(fig[1, 1],title="μ")
    ax2 = Axis(fig[1, 2], title = "π")
    ax3 = Axis(fig[2, 1], title = "Cross Sectional μ")
    ax4 = Axis(fig[2, 2], title = "Cross Sectional π")

    # Create lines for each plot
    firm_density_line = lines!(ax, p.x, sol(0), linewidth = 3, label = "Firms' density")
    profits_line = lines!(ax2, p.x, zeros(length(p.x)), linewidth = 3, label = "Profits per sector")
    firms_crossectional = lines!(ax3, crossμX0, crossμY, linewidth = 3, label = "Cross sectional firms' density")
    profits_crossectional = lines!(ax4, crossProfitsX0, crossProfitsY, linewidth = 3, label = "Cross sectional profits")

    ax.limits = (0, n, 0.9*minimum(sol(0)), 1.1*maximum(sol(0)))
    ax2.limits = (0, n, 0.9*minimum(profits), 1.1*maximum(profits))
    ax3.limits = (minimum(crossμX0), maximum(crossμX0), 0.9*minimum(crossμY), 1.1*maximum(crossμY))
    ax4.limits = (minimum(filter!(!isnan,crossProfitsX0)), maximum(filter!(!isnan,crossProfitsX0)), 0.9*minimum(filter!(!isnan,crossProfitsY)), 1.1*maximum(filter!(!isnan,crossProfitsY)))
    display(fig)

    record(fig, fName, LinRange(0, T_end, NSteps)) do t
        @show t

        profits = computeProfits(sol(t),p)
        profits[isinf.(profits)] .= 0    
        crossμX,crossμY = toCrossectional(sol(t),p)
        crossProfitsX,crossProfitsY = toCrossectional(profits,p)
    
        # Update the data in the lines
        firm_density_line[2] = sol(t)
        profits_line[2] = profits
        firms_crossectional[1] = crossμX
        firms_crossectional[2] = crossμY
        profits_crossectional[1] = crossProfitsX
        profits_crossectional[2] = crossProfitsY

        # production_line[2] = production
        # derivative_profits_line[2] = ∂x(profits,p)
        # smoothed_derivative_profits_line[2] = convolve(∂x(profits,p),K,p)

        # Update title
        ax.title = "μ, t = " * string(round(t, digits = 3))
        ax2.title = "π, t = " * string(round(t, digits = 3))
        ax3.title = "Cross Sectional μ, t = " * string(round(t, digits = 3))
        ax4.title = "Cross Sectional π, t = " * string(round(t, digits = 3))

        ax.limits = (0, n, 0.9*minimum(sol(t)), 1.1*maximum(sol(t)))
        ax2.limits = (0, n, 0.9*minimum(profits), 1.1*maximum(profits))
        # ax3.limits = (minimum(crossμX0), maximum(crossμX0), 0.9*minimum(crossμY), 1.1*maximum(crossμY))
        # ax4.limits = (minimum(crossProfitsX0), maximum(crossProfitsX0), 0.9*minimum(crossProfitsY), 1.1*maximum(crossProfitsY))

        display(fig)
        sleep(0.01)
    end

end

"""
profitsTwoTimes(sol, p)

- Plots cross-sectional profits at t=0 and t=T_end.
"""
function profitsTwoTimes(sol,p)
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1],title="π")

    crossProfitsX0,crossProfitsY0 = toCrossectional(computeProfits(sol(0),p),p)
    crossProfitsXT,crossProfitsYT = toCrossectional(computeProfits(sol(p.T_end),p),p)

    lines!(ax, crossProfitsX0, crossProfitsY0, linewidth = 3, label = "Profits at period 0")
    lines!(ax, crossProfitsXT, crossProfitsYT, linewidth = 3, label = "Profits at period T")

    axislegend(ax)
    display(fig)
end

function profitsVSΔProfits(sol,p)
    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1])

    profits0  = computeProfits(sol(0), p)
    profitsT = computeProfits(sol(p.T_end), p)
    Δprofits = profitsT - profits0

    plot!(ax, profits0, Δprofits,label="")
    ax.limits = (0.9*minimum(profits0), 1.1*maximum(profits0), 0.9*minimum(Δprofits), 1.1*maximum(Δprofits))
    axislegend(ax)
    display(fig)

end

function profitsVSΔProfitsOverTime(sol,p; NSteps = 100,fName = "profitsVSDeltaProfits.mp4")

    fig = Figure(size = (1000, 1000))
    ax = Axis(fig[1, 1], xlabel = "π(0)", ylabel = "π(t)-π(0)", title = "")

    profits0  = computeProfits(sol(0), p)
    profitsT = computeProfits(sol(p.T_end), p)
    Δprofits = profitsT - profits0
    ax.limits = (0.9*minimum(profits0), 1.1*maximum(profits0), 0.9*minimum(Δprofits), 1.1*maximum(Δprofits))

    scatterProfits = scatter!(ax, profits0, zeros(length(profits0)),label="")
    display(fig)

    record(fig, fName, LinRange(0, p.T_end, NSteps)) do t
        
        @show t
        profitst = computeProfits(sol(t), p)
        Δprofits = profitst - profits0
        scatterProfits[2] = Δprofits
        
        # Compute profits
        ax.title = "t = " * string(round(t, digits = 3))
        axislegend(ax)
        
        display(fig)
        sleep(0.01)
    end


end


