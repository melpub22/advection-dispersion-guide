% =========================================================
% Advection Dispersion Equations Implemented in 1 Dimension
% ---------------------------------------------------------
% Author:   github.com/melpub22
% Creation: 26 November 2023
% License:  MIT License
% ---------------------------------------------------------
% Implementation of a simple fixed-bed advection-dispersion
% scenario in one dimension with MATLAB's pdepe function.
%
% SCENARIO: Consider a fixed bed of length L = 10 m with a
% constant inlet concentration C_in = 100 mg/L. Initially,
% the bed contains no fluid. Fluid begins to flow into the
% bed with a velocity of v = 1e-3 m/s and the dispersion
% coefficient is D = 2e-2 m^2/s.
%
% OUTPUTS: A plot of concentration profiles in the bed over
% time and a breakthrough curve.
%
% LIMITATIONS: This is a single-species scenario, and only
% advection and dispersion are considered. The boundary
% conditions set for this scenario may not be applicable to
% your scenario.
% =========================================================

% NOTE: This same scenario is implemented in Julia at:
% https://github.com/melpub22/advection-dispersion-guide/tree/main/AdvectionDispersion1D/julia_ade_1d.jl

t_max = 2000;   % Duration of study, s
L_max = 10;     % Length of fixed bed, m
 
t_num = 500;    % Number of values in time grid
L_num = 100;    % Number of values in position grid
 
% Create linear grids for time and position
t = linspace(0, t_max, t_num);
x = linspace(0, L_max, L_num);
 
% Solve the PDE system
sol = pdepe(0, @pde_cfs, @pde_ic, @pde_bc, x, t);
 
% Plot the concentration profile at various times
figure;
plot(x, sol(100, :))
hold on
plot(x, sol(200, :))
plot(x, sol(300, :))
plot(x, sol(400, :))
plot(x, sol(500, :))
title("Concentration profiles over time")
xlabel("Position x, m")
ylabel("Concentration u_x, mg/L")
legend("t = 400 s", "t = 800 s", "t = 1200 s", "t = 1600 s", "t = 2000 s", "Location", "northeast")
axis([0 L_max, 0 100])
 
% Plot the breakthrough curve
figure;
plot(t, sol(:, end))
title("Breakthrough curve")
xlabel("Time t, s")
ylabel("Outlet concentration u_r, mg/L")
 
function [c, f, s] = pde_cfs(x, t, u, DuDx)
    D = 2e-2;           % Dispersion coefficient, m^2/s
    v = 1e-3;           % Velocity, m/s
 
    c = 1;              % c-function required for pdepe()
    f = D*DuDx - v*u;   % f-function required for pdepe()
    s = 0;              % s-function required for pdepe()
end
 
function ic = pde_ic(x)
    ic = 0;             % Initial concentration at each position, mg/L
end
 
function [pl, ql, pr, qr] = pde_bc(xl, ul, xr, ur, t)
    v = 1e-3;           % Velocity, m/s
    Cin = 100;          % Concentration at inlet, mg/L
 
    pl = ul - Cin;      % p-function for pdepe() boundary condition on the left (inlet)
    ql = 0;             % q-function for pdepe() boundary condition on the left (inlet)
    pr = v*ur;          % p-function for pdepe() boundary condition on the right (outlet)
    qr = 1;             % q-function for pdepe() boundary condition on the right (outlet)
end
