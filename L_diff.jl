

function df!(du,u,p,t)
    @unpack_parameters p
    Integral = sum(Aβρ .* u.^((1+η/(1-β))*(ρ-1)/ρ)) * Δx
    du .= -(1-β)/Integral * ∂x(u .* ∂x(Aβρ .* u.^(η/(1-β)*(ρ-1)/ρ - 1/ρ),p),p)       
    @show t
end

function dfFixedCost!(du,u,p,t)
    @unpack_parameters p
    Integral = sum(Aβρ .* u.^((1+η/(1-β))*(ρ-1)/ρ)) * Δx
    ∂xπ = ∂x(Aβρ .* u.^(η/(1-β)*(ρ-1)/ρ - 1/ρ),p)

    du .= -(1-β)/Integral * ∂x(u .* ∂xπ,p) .* convolve(Float64.(abs.(∂xπ) .> c₀),K,p)
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


