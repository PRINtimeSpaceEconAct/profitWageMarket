function computeProfits(μ,p)
    @unpack_parameters p
    profits = (1-β) * ALβσ .* μ.^((η-β)*(σ-1)/σ - 1/σ) / sum(ALβσ .* μ.^((1+η-β)*(σ-1)/σ) * Δx)
    return profits
end

function computeProfitsMobile(μ,p)
    @unpack_parameters p
    α = (σ - 1) / (σ*(1 - β) + β)
    γ = (η*(σ - 1) - 1) / (σ*(1 - β) + β)
    profits = (1 - β) * A.^α .* μ.^γ / sum((A .* μ.^(1 - β + η)).^α * Δx)
    return profits
end

function computeProduction(μ,p)
    @unpack_parameters p
    Y = ALβσ .* μ.^((η-β)*(σ-1)/σ) / sum(ALβσ .* μ.^((1+η-β)*(σ-1)/σ) * Δx)
    return Y
end

function computeP(μ,p)
    @unpack_parameters p
    P = sum((A .* L.^β .* μ.^(1+η-β)).^(-(1-σ)/σ) * Δx).^(-σ/(σ-1))
    return P
end

function computePMobile(μ,p)
    @unpack_parameters p
    P = sum((A .* μ.^(1-β)).^((σ-1)/(σ*(1-β)+β)) * Δx)^((σ*(1-β)+β)/(1-σ))
    return P
end


function computeX(μ,p) 
    X = 1 ./ computeP(μ,p)
    return X
end

function computeXMobile(μ,p) 
    @unpack_parameters p
    α = (σ - 1) / (σ*(1 - β) + β)
    X = sum((A .* μ.^(1 - β + η)).^α * Δx) ^ ((σ*(1 - β) + β) / (σ - 1))
    return X
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



# function doIt(f,p)

#     # f(x) = -tan(5.5*(x-1/4))+30
    
#     y = f.(p.x[1:round(Int,p.Nx/2)])
#     μ = mirror(y)
#     μ = μ/sum(μ*p.Δx)
#     # μ = convolve(μ,K,p)
#     profits = computeProfits(μ,p)
#     profits[isinf.(profits)] .= 0    
#     crossProfitsX,crossProfitsY = toCrossectional(profits,p)
#     crossμX,crossμY = toCrossectional(μ,p)
#     crossμXDiscrete,crossμYDiscrete = toCrossectional(μ)
#     crossProfitsXDiscrete,crossProfitsYDiscrete = toCrossectional(profits)


#     fig = Figure(size=(1000,1000))
#     ax = Axis(fig[1, 1],title="μ")
#     lines!(ax,p.x,μ)

#     ax2 = Axis(fig[1, 2], title = "π")
#     lines!(ax2,p.x,profits)

#     ax3 = Axis(fig[2, 1], title = "Cross Sectional μ")
#     lines!(ax3,crossμX,crossμY)
#     # lines!(ax3,crossμXDiscrete,crossμYDiscrete)

#     ax4 = Axis(fig[2, 2], title = "Cross Sectional π")
#     lines!(ax4,crossProfitsX,crossProfitsY)
#     # lines!(ax4,crossProfitsXDiscrete,crossProfitsYDiscrete)

#     display(fig)

#     #save figure
#     save("CrossSectionalProfits.png",fig)
# end

