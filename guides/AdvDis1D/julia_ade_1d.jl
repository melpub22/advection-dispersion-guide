# =========================================================
# Advection Dispersion Equations Implemented in 1 Dimension
# ---------------------------------------------------------
# Author:   github.com/melpub22
# Creation: 26 November 2023
# License:  MIT License
# ---------------------------------------------------------
# Implementation of a simple fixed-bed advection-dispersion
# scenario in one dimension with Julia.
#
# SCENARIO: Consider a fixed bed of length L = 10 m with a
# constant inlet concentration C_in = 100 mg/L. Initially,
# the bed contains no fluid. Fluid begins to flow into the
# bed with a velocity of v = 1e-3 m/s and the dispersion
# coefficient is D = 2e-2 m^2/s.
#
# OUTPUTS: A plot of concentration profiles in the bed over
# time and a breakthrough curve.
#
# LIMITATIONS: This is a single-species scenario, and only
# advection and dispersion are considered. The boundary
# conditions set for this scenario may not be applicable to
# your scenario.
# =========================================================

# NOTE: This same scenario is implemented in MATLAB at:
# https://github.com/melpub22/advection-dispersion-guide/tree/main/AdvectionDispersion1D/matlab_ade_1d_pdepe.m

using DomainSets, ModelingToolkit, MethodOfLines, OrdinaryDiffEq, Plots

# Define scenario-specific variables
D = 2.0e-2      # Dispersion coefficient, m²/s
v = 1.0e-3      # Velocity, m/s
L = 10.0        # Bed length, m
t_max = 2000.0  # Simulation duration, s
Cin = 100.0     # Inlet concentration, mg/L

t_num = 500     # Number of values in the time grid
L_num = 100     # Number of values in the position grid

# Setup parameters and variables for the PDE
@parameters t x
@variables u(..)
Dt = Differential(t)
Dx = Differential(x)
Dxx = Differential(x)^2

# Create the Advection-Dispersion PDE
pde = Dt(u(t, x)) ~ D*Dxx(u(t, x)) - v*Dx(u(t, x))

# Define the initial and boundary conditions
bcs =  [u(0, x) ~ 0,            # Initial condition requiring that u = 0 at all x for t = 0
        u(t, 0) ~ Cin,          # Boundary condition holding inlet concentration constant with u = Cin at x = 0
        Dx(u(t, L)) ~ 0]        # Boundary condition requiring that ∂u/∂x = 0 at x = L

# Define time and position domains
domains =  [t ∈ Interval(0.0, t_max),   # Time domain
            x ∈ Interval(0.0, L)]       # Position domain

# Create the PDE system
@named pde_system = PDESystem(pde, bcs, domains, [t, x], [u(t, x)])

# Define the time and position deltas
dt = t_max / t_num              # Time delta
dx = L / L_num                  # Position delta

# Create the methods of lines discretization
discretization = MOLFiniteDifference([x => dx], t)

# Setup an ODE problem from the methods of lines descretization
problem = discretize(pde_system, discretization)

# Solve the PDE system
solution = solve(problem, Tsit5(), saveat=dt)

# Extract the discrete x and t values and the solution u (concentration) values
discrete_x = solution[x]        # Extract the discrete x values
discrete_t = solution[t]        # Extract the discrete t values
solution_u = solution[u(t, x)]  # Extract the solution u values

# Plot the profiles curve
profiles = plot()
plot!(discrete_x, solution_u[101, :], label="t = 400 s")
plot!(discrete_x, solution_u[201, :], label="t = 800 s")
plot!(discrete_x, solution_u[301, :], label="t = 1200 s")
plot!(discrete_x, solution_u[401, :], label="t = 1600 s")
plot!(discrete_x, solution_u[501, :], label="t = 2000 s")
title!("Concentration profiles over time")
xlabel!("Position x, m")
ylabel!("Concentration uₓ, mg/L")
display(profiles)

# Plot the breakthrough curve
breakthrough = plot(discrete_t, solution_u[:, end], label="Breakthrough concentration")
title!("Breakthrough curve")
xlabel!("Time t, s")
ylabel!("Outlet concentration uᵣ, mg/L")
display(breakthrough)

# This script follows the methods of lines approach outlined by the 
# MethodOfLines.jl heat equation tutorial available at:
# https://github.com/SciML/MethodOfLines.jl/blob/master/docs/src/tutorials/heat.md
# MethodOfLines.jl and its tutorials has the following license:
#        MIT License
#
#        Copyright (c) 2022 SciML Open Source Scientific Machine Learning Organization
#
#        Permission is hereby granted, free of charge, to any person obtaining a copy
#        of this software and associated documentation files (the "Software"), to deal
#        in the Software without restriction, including without limitation the rights
#        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#        copies of the Software, and to permit persons to whom the Software is
#        furnished to do so, subject to the following conditions:
#
#        The above copyright notice and this permission notice shall be included in all
#        copies or substantial portions of the Software.
#
#        THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#        SOFTWARE.
