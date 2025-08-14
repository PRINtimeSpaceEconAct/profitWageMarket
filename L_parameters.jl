
@with_kw mutable struct parameters

    # model     
    β::Float64 = 0.86
    # β::Float64 = 0.5
    ρ1::Float64 = 0.214
    ρ2::Float64 = 0.04
    η::Float64 = β * (1+ρ2/ρ1) - 1  
    # η::Float64 = 1.2
    σ::Float64 = 1 + ρ1 / (β*(1+ρ2) - ρ1)
    # σ::Float64 = 2.0
    ρ::Float64 = σ*(1-β) + β
    # c₀::Float64 = 0.05
    c₀::Float64 = 0.1



    # domain
    n::Float64 = 1.0
    # T_end::Float64 = 2.0
    T_end::Float64 = 2.0
    
    # numerical
    Δx::Float64 = 5*1e-3
    Nx::Int = Int(n/Δx)
    x::LinRange{Float64, Int64} = LinRange(0,n-Δx,Nx)
    t_save::LinRange{Float64, Int64} = LinRange(0,T_end,1000)
    h::Float64 = 0.05    # smoothing parameter when needed

    # initial condition
    freq::Int64 = 1
    height::Float64 = 1.5
    mass::Float64 = 1.0
    shift::Float64 = 0.25
    μ₀::Vector{Float64} = initialCondition(Δx,x,n,freq,height,mass,shift)

    # technological progress
    shiftA::Float64 = 0.0
    # A::Vector{Float64} = computeA(Δx,x,n,4,height,2.0,-0.25)
    A::Vector{Float64} = ones(Nx)
    # L::Vector{Float64} = ones(Nx)
    L::Vector{Float64} = computeA(Δx,x,n,1,height,1.0,-0.25)
    ALβσ::Vector{Float64} = (A.*L.^β).^((σ-1)/σ)

    # checks
    @assert σ != 1.0
    @assert height > 1
end

function initialCondition(Δx,x,n,freq,height,mass,shift)
    """ initialize initial condition for PDE by sinus"""

    Nx = length(x)
    f(x) = -tan(5.5*(x-1/4))+30
    y = f.(x[1:round(Int,Nx/2)])
    μ₀ = mirror(y)
    μ₀ = μ₀/sum(μ₀*Δx)

    # μ₀ .= shuffle(μ₀)
    # circshift!(μ₀,round(Int,Nx/5))
    # μ₀ = sin.(2*π*(x .+ shift)/n*freq ) .+ height
    # μ₀ = mass * μ₀ / ( sum(μ₀) * Δx )
    return μ₀
end

function computeA(Δx,x,n,freq,height,mass,shiftA)
    """ define the technological progress A(i)"""
    A = 0.5*sin.(2*π*(x .+ shiftA)/n*freq) .+ height
    A = mass * A / ( sum(A) * Δx )
    return A
end

