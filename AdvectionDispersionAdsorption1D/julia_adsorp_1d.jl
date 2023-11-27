# =========================================================
# Adsorption column PDEs implemented in 1 Dimension
# ---------------------------------------------------------
# Author:   github.com/melpub22
# Creation: 27 November 2023
# License:  MIT License
# ---------------------------------------------------------
# Implementation of a simple fixed-bed advection-dispersion
# scenario in one dimension with Julia.
#
# SCENARIO: Consider a fixed bed of length L = 10 m with a
# constant inlet concentration C_in = 100 mol/m^3.
# Initially, the bed contains no fluid. Fluid begins to
# flow into the bed with a velocity of v = 1e-3 m/s and the
# dispersion coefficient is D = 2e-2 m^2/s. 
#
# The fixed bed has a porosity of 0.3 and the Langmuir
# adsorption isotherm model parameters: max sorbent loading
# capacity q_max = 1.5e-5 mol/kg, Langmuir coefficient
# K_L = 1e3 m^3/mol, and density of sorbent in column
# rho = 1e6 kg/m^3.
#
# The linear driving force coefficient K_F = 0.1 s^-1.
#
# OUTPUTS: A plot of concentration profiles in the bed over
# time and a breakthrough curve.
#
# LIMITATIONS: This is a single-species scenario, and only
# advection, dispersion, and adsorptinn are considered. The
# boundary conditions set for this scenario may not be
# applicable to your scenario.
# =========================================================

# Learn more about the MethodOfLines.jl package that this
# script depends on at: docs.sciml.ai/MethodOfLines/stable/

using DomainSets, ModelingToolkit, MethodOfLines, OrdinaryDiffEq, Plots

# Define advection-dispersion model parameters
D = 2.0e-2      # Dispersion coefficient, m²/s
v = 1.0e-3      # Velocity, m/s
L = 10.0        # Bed length, m
Cin = 100.0     # Inlet concentration, mol/m³

# Define additional parameters for adsorption
ε = 0.3         # Porosity
ρ = 1.0e6       # Density, kg/m³
K_F = 0.1       # Linear driving force coefficient, s⁻¹

# Define adsorption isotherm parameters
q_max = 1.5e-5  # Maximum sorbent loading, mol/kg
K_L = 1.0e3     # Langmuir coefficient, m³/mol

# Other simulation parameters
t_max = 2000.0  # Simulation duration, s
t_num = 500     # Number of values in the time grid
L_num = 100     # Number of values in the position grid

# Define the Langmuir adsorption isotherm model function
q_star(c) = (q_max * K_L * c)/(1 + K_L * c)

# Setup parameters and variables for the PDE
@parameters t x
@variables C(..) q(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# Create the Advection-Dispersion PDE
pde =  [Dt(C(t, x)) ~ D*Dxx(C(t, x)) - v*Dx(C(t, x)) - ((1-ε)/ε)*ρ*Dt(q(t, x)),
        Dt(q(t, x)) ~ K_F*(q_star(C(t, x)) - q(t,x))]

# Define the initial and boundary conditions
bcs =  [C(0, x) ~ 0,            # Initial condition requiring that C = 0 at all x for t = 0
        C(t, 0) ~ Cin,          # Boundary condition holding inlet concentration constant with C = Cin at x = 0
        Dx(C(t, L)) ~ 0,        # Boundary condition requiring that ∂C/∂x = 0 at x = L
        q(0, x) ~ 0,            # Initial condition requiring no sorbent loading at t = 0
        Dx(q(t, 0)) ~ 0,        # Boundary condition requiring that ∂q/∂x = 0 at x = 0
        Dx(q(t, L)) ~ 0]        # Boundary condition requiring that ∂q/∂x = 0 at x = L

# Define time and position domains
domains =  [t ∈ Interval(0.0, t_max),   # Time domain
            x ∈ Interval(0.0, L)]       # Position domain

# Create the PDE system
@named pde_system = PDESystem(pde, bcs, domains, [t, x], [C(t, x), q(t, x)])

# Define the time and position deltas
dt = t_max / t_num              # Time delta
dx = L / L_num                  # Position delta

# Create the methods of lines discretization
discretization = MOLFiniteDifference([x => dx], t)

# Setup an ODE problem from the methods of lines descretization
problem = discretize(pde_system, discretization)

# Solve the PDE system
solution = solve(problem, QNDF(), saveat=dt)

# Extract the discrete x and t values and the solution values
discrete_x = solution[x]        # Extract the discrete x values
discrete_t = solution[t]        # Extract the discrete t values
solution_C = solution[C(t, x)]  # Extract the solution C values
solution_q = solution[q(t, x)]  # Extract the solution q values

# Plot the profiles curve
profiles = plot()
plot!(discrete_x, solution_C[101, :], label="t = 400 s")
plot!(discrete_x, solution_C[201, :], label="t = 800 s")
plot!(discrete_x, solution_C[301, :], label="t = 1200 s")
plot!(discrete_x, solution_C[401, :], label="t = 1600 s")
plot!(discrete_x, solution_C[501, :], label="t = 2000 s")
title!("Concentration profiles over time")
xlabel!("Position x, m")
ylabel!("Concentration uₓ, mol/m³")
display(profiles)

# Plot the breakthrough curve
breakthrough = plot(discrete_t, solution_C[:, end], label="Breakthrough concentration")
title!("Breakthrough curve")
xlabel!("Time t, s")
ylabel!("Outlet concentration uᵣ, mol/m³")
display(breakthrough)

# Plot the sorbent loading
loading = plot()
plot!(discrete_x, solution_q[26, :], label="t = 100 s")
plot!(discrete_x, solution_q[51, :], label="t = 200 s")
plot!(discrete_x, solution_q[101, :], label="t = 400 s")
plot!(discrete_x, solution_q[201, :], label="t = 800 s")
plot!(discrete_x, solution_q[251, :], label="t = 1000 s")
title!("Sorbent loading over time")
xlabel!("Position x, m")
ylabel!("Sorbent loading q, mol/kg")
display(loading)