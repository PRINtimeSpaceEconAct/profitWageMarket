function computeProfits(μ,p)
    @unpack_parameters p
    profits = (1-β) * Aβρ .* μ.^(η/(1-β)*(ρ-1)/ρ - 1/ρ) / sum(Aβρ .* μ.^((1+η/(1-β))*(ρ-1)/ρ) * Δx)
    return profits
end

function computeLabour(μ,p)
    @unpack_parameters p
    L = (A .* μ.^(1-β+η)).^((σ-1)/ρ) / sum(Aβρ .* μ.^( (1+η/(1-β)) * (ρ-1)/ρ) * Δx)
    return L
end

function computeProduction(μ,p)
    Y = computeLabour(μ,p)
    return Y
end

function toCrossectional(f;Np = 100)
    y = LinRange(minimum(f),maximum(f),Np)
    yCrossectional = zeros(Np-1)

    for i in 1:(Np-1)
        yCrossectional[i] = sum( (f .> y[i]) .&  (f .<= y[i+1]) ) / length(f)
    end

    return y[2:end],yCrossectional
end


function toCrossectional(f,p)
    y = similar(f)
    ∂xf = ∂x(f,p)
    for (i,fi) in enumerate(f)
        y[i] = abs(1/∂xf[i])
    end

    p = sortperm(f)
    fPerm = f[p]
    yPerm = y[p]

    # normalize
    yPerm .= yPerm / sum( yPerm[2:end] .* diff(fPerm) )
    return fPerm,yPerm
end


function K(x)
    """ smoothing kernel """
    return abs(x) <= 1 ? 1-abs(x) : 0.0
    # return 1/√π * exp(-x^2)
end

function K(x,h)
    """ rescaled kernel """
    return 1/h*K(x/h)
end

# distance on the 1D torus:
function dist(x,y,L)
    return min( abs(x-y), L - abs(x-y) )
end

#one dimensional convolution on the thorus:
function convolve(f,K,p)
    """ convolution of f and g on the torus"""
    N = length(f)
    fConv = similar(f)
    for i in 1:N
        fConv[i] = sum( [f[j] * K(dist(p.x[i],p.x[j],p.n),p.h) * p.Δx for j in 1:N] )
    end
    return fConv
end

function mirror(f)
    """ mirror the function f"""
    return [f;f[end:-1:1]]
end

function doIt(f,p)

    # f(x) = -tan(5.5*(x-1/4))+30
    
    y = f.(p.x[1:round(Int,p.Nx/2)])
    μ = mirror(y)
    μ = μ/sum(μ*p.Δx)
    # μ = convolve(μ,K,p)
    profits = computeProfits(μ,p)
    profits[isinf.(profits)] .= 0    
    crossProfitsX,crossProfitsY = toCrossectional(profits,p)
    crossμX,crossμY = toCrossectional(μ,p)
    crossμXDiscrete,crossμYDiscrete = toCrossectional(μ)
    crossProfitsXDiscrete,crossProfitsYDiscrete = toCrossectional(profits)


    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1, 1],title="μ")
    lines!(ax,p.x,μ)

    ax2 = Axis(fig[1, 2], title = "π")
    lines!(ax2,p.x,profits)

    ax3 = Axis(fig[2, 1], title = "Cross Sectional μ")
    lines!(ax3,crossμX,crossμY)
    # lines!(ax3,crossμXDiscrete,crossμYDiscrete)

    ax4 = Axis(fig[2, 2], title = "Cross Sectional π")
    lines!(ax4,crossProfitsX,crossProfitsY)
    # lines!(ax4,crossProfitsXDiscrete,crossProfitsYDiscrete)

    display(fig)

    #save figure
    save("CrossSectionalProfits.png",fig)
end

