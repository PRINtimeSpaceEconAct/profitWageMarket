

function df!(du,u,p,t)
    @unpack_parameters p
    profits = computeProfits(u,p)
    du .= -(1-β) * ∂x(u .* ∂x(profits,p),p)       
    @show t
end

function dfFixedCost!(du,u,p,t)
    @unpack_parameters p
    profits = computeProfits(u,p)
    ∂xπ = ∂x(profits,p)
    
    du .= -(1-β) * ∂x(u .* ∂xπ,p) .* convolve(Float64.(abs.(∂xπ) .> c₀),K,p)
    # du .= -(1-β) * ∂x(u .* ∂xπ,p) .* Float64.(abs.(∂xπ) .> c₀)
    @show t
end


function Δ(f,p)
    return Δ(f,p.Nx,p.Δx)
end

function Δ(f,Nx,Δx)
    """ laplacian periodic boundary"""
    Δf = similar(f)
    for i=2:Nx-1
        Δf[i] = (f[i+1] - 2*f[i] + f[i-1])/Δx^2
    end
    Δf[1] = (f[2] - 2*f[1] + f[Nx])/Δx^2
    Δf[Nx] = (f[1] - 2*f[Nx] + f[Nx-1])/Δx^2
    return Δf
end

function ∂x(f,p)
    return ∂x(f,p.Nx,p.Δx)
end

function ∂x(f,Nx,Δx)
    """ derivative periodic boundary"""
    ∂f = similar(f)
    for i=2:Nx-1
        ∂f[i] = (f[i+1] - f[i-1])/(2*Δx)
    end
    ∂f[1] = (f[2] - f[Nx])/(2*Δx)
    ∂f[Nx] = (f[1] - f[Nx-1])/(2*Δx)
    return ∂f
end


