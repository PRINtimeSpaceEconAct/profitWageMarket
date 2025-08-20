# Code used in ''The invisible hand as an emergent property: a gradient flow approach''

A Julia codebase for simulating firm density dynamics over a continuous goods space with heterogeneous technology and labor.

## Quick Start

```julia
using Pkg
Pkg.activate(@__DIR__)

# Run all models and generate figures
p = parameters()
runAll(p)
```

## Core Models

- **Baseline**: Standard profit-driven firm entry/exit dynamics
- **Mobile**: Alternative profit and price index formulation  
- **Fixed Cost**: Includes fixed costs in firm dynamics

## Key Files

- `main.jl` - Main solver functions and `runAll()` convenience function
- `L_parameters.jl` - Model parameters and grid setup
- `L_utils.jl` - Core economic functions (profits, prices, aggregates)
- `L_plots_paper.jl` - Publication-ready figure generation
- `L_plots.jl` - Interactive plotting and animations

## Model Mechanics

The state variable is firm density $\mu(x,t)$ over goods space $x \in [0,n]$. At each point:

- Technology levels A(x) and labor L(x) determine local productivity
- Firms reallocation based on profit differentials
- Aggregate price index P depends on the full distribution
- Cross-sectional analysis transforms spatial to ranked distributions

## Usage Examples

```julia
# Solve individual models
p = parameters(n = 1.0, T_end = 2.0)
sol, p = solveModel(p)          # Baseline
sol_mob, p = solveModelMobile(p) # Mobile variant  
sol_fc, p = solveModelFixedCost(p) # Fixed cost

# Generate specific figures
plot1(sol, p, 0.0, 0.2)  # Density and profits over space
plot2(sol, p, 0.01, 0.5, 2.0)  # Profit variation scatter
# others related to other figures
```

## Customization

- Modify parameters in `L_parameters.jl`

