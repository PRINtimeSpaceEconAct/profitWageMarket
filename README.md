# Profit–Wage Market Dynamics (Julia)

## Overview
Codebase to simulate spatial/continuous–variety firm density dynamics with heterogeneous technology and labor, and to produce paper–ready figures. The state variable is the firms' density μ(x,t) over a 1D goods space x ∈ [0,n]. Profits, production, price index and aggregate variables are computed from μ. Several model variants:
- Baseline (free entry diffusion / selection)
- Mobile variant (alternative profit / price index definitions)
- Fixed cost variant

## Key Mechanisms (selected formulas)
Let β, σ, η be parameters and A(x), L(x) profiles.
- Profits (baseline):
  π(i) = (1-β) * (A(i) L(i)^β)^{(σ-1)/σ} * μ(i)^{((η-β)(σ-1)/σ - 1/σ)} / ∫ (A L^β)^{(σ-1)/σ} μ^{(1+η-β)(σ-1)/σ} dx
- Profits (mobile):
  π(i) = (1-β) * A(i)^α μ(i)^γ / ∫ [A(j) μ(j)^{1-β+η}]^α dj, with α = (σ-1)/(σ(1-β)+β), γ = (η(σ-1)-1)/(σ(1-β)+β)
- Price index (baseline):
  P = [∫ (A L^β μ^{1+η-β})^{-(1-σ)/σ} dx]^{-σ/(σ-1)}
- Price index (mobile):
  P = [∫ (A μ^{1-β})^α dx]^{-(σ(1-β)+β)/(σ-1)}
- Aggregate X = 1/P
- Cross–sectional distributions built from sorted values and inverse gradient density weighting.

(See L_utils.jl for precise implementations.)

## File Structure
- main.jl: Entry points solveModel, solveModelMobile, solveModelFixedCost (sets up and solves ODE problems).
- L_LoadAll.jl: Central include/use (loads parameters, utilities, plots, ODE definitions).
- L_parameters.jl: Parameter struct (domain, numerical grid, initial condition, technology profiles).
- L_utils.jl: Core economic functions (profits, production, price indices, aggregates, cross-sectional transforms, kernels).
- L_plots.jl: Exploratory dynamic and cross–sectional plotting (animations).
- L_plots_paper.jl: Publication–style figures (monochrome theme, line styles).
- (Implied) ODE definition files (e.g. df!, dfMobile!, dfFixedCost!) define μ dynamics.
- Project.toml / Manifest.toml (if present): environment specification.

## Typical Workflow
1. Activate environment:
   using Pkg; Pkg.activate(@__DIR__)
2. Set / modify parameters:
   p = parameters(β=0.86)
3. Solve baseline:
   sol,p = solveModel(p)
4. Generate paper figures:
   include("L_plots_paper.jl"); plot1(sol,p,0.0,0.2)
5. Solve alternative variants:
   solMob,p = solveModelMobile(p)
   solFC,p = solveModelFixedCost(p)
6. Inspect diagnostics or export animations via plot functions in L_plots.jl.

## Numerical Notes
- Space discretized by uniform grid x with step Δx.
- Time solved via DifferentialEquations.jl (Rosenbrock23, semi-implicit/stiff-friendly).
- Integrals approximated by Riemann sums (sum( f .* Δx )).
- Cross-sectional transformations rely on sorting and derivative-based density weighting.

## Extending
- Add new profit specification: implement in L_utils.jl and mirror pattern of computeProfits*.
- Add new model dynamics: define dfNew! and a corresponding solveModelNew wrapper.
- Add figure: follow style in L_plots_paper.jl (monochrome, line styles, axislegend).

## Reproducibility
Set random seeds (if randomness added) before initial condition. Version pin packages via Project/Manifest.

## Citation
If used in academic work, cite relevant theoretical paper plus this repository.

