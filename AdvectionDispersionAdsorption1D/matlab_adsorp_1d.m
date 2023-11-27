% =========================================================
% Adsorption column PDEs implemented in 1 Dimension
% ---------------------------------------------------------
% Author:   github.com/melpub22
% Creation: 27 November 2023
% License:  MIT License
% ---------------------------------------------------------
% Implementation of a simple fixed-bed advection-dispersion
% scenario in one dimension with MATLAB's pdepe function.
%
% SCENARIO: Consider a fixed bed of length L = 10 m with a
% constant inlet concentration C_in = 100 mol/m^3.
% Initially, the bed contains no fluid. Fluid begins to
% flow into the bed with a velocity of v = 1e-3 m/s and the
% dispersion coefficient is D = 2e-2 m^2/s. 
%
% The fixed bed has a porosity of 0.3 and the Langmuir
% adsorption isotherm model parameters: max sorbent loading
% capacity q_max = 1.5e-5 mol/kg, Langmuir coefficient
% K_L = 1e3 m^3/mol, and density of sorbent in column
% rho = 1e6 kg/m^3.
%
% The linear driving force coefficient K_F = 0.1 s^-1.
%
% OUTPUTS: A plot of concentration profiles in the bed over
% time and a breakthrough curve.
%
% LIMITATIONS: This is a single-species scenario, and only
% advection, dispersion, and adsorptinn are considered. The
% boundary conditions set for this scenario may not be
% applicable to your scenario.
% =========================================================

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
plot(x, sol(100, :, 1))
hold on
plot(x, sol(200, :, 1))
plot(x, sol(300, :, 1))
plot(x, sol(400, :, 1))
plot(x, sol(500, :, 1))
title("Concentration profiles over time")
xlabel("Position x, m")
ylabel("Concentration u_x, mol/m^3")
legend("t = 400 s", "t = 800 s", "t = 1200 s", "t = 1600 s", "t = 2000 s", "Location", "northeast")
axis([0 L_max, 0 100])
 
% Plot the breakthrough curve
figure;
plot(t, sol(:, end, 1))
title("Breakthrough curve")
xlabel("Time t, s")
ylabel("Outlet concentration u_r, mol/m^3")

% Plot the adsorbent loading curve
figure;
plot(x, sol(25, :, 2))
hold on
plot(x, sol(50, :, 2))
plot(x, sol(100, :, 2))
plot(x, sol(200, :, 2))
plot(x, sol(250, :, 2))
title("Sorbent loading over time")
xlabel("Position x, m")
ylabel("Sorbent loading q, mol/kg")
legend("t = 100 s", "t = 200 s", "t = 400 s", "t = 800 s", "t = 1000 s", "Location", "northeast")

function [c, f, s] = pde_cfs(x, t, u, DuDx)
    D = 2e-2;                       % Dispersion coefficient, m^2/s
    v = 1e-3;                       % Velocity, m/s
    por = 0.3;                      % Porosity
    rho = 1e6;                      % Density, kg/m^3
    K_F = 0.1;                      % Linear driving force coefficient, s^-1

    % Relationship for the linear driving force model
    DqDt = K_F * (isotherm(u(1)) - u(2));
 
    % Handling the first PDE: The Advection-Diffusion equation
    c(1) = 1;                       % c-function required for pdepe()
    f(1) = D*DuDx(1) - v*u(1);      % f-function required for pdepe()
    s(1) = -((1-por)/por)*rho*DqDt; % s-function required for pdepe()

    % Handling the second PDE: The Linear Driving Force equation
    c(2) = 1;                       % c-function required for pdepe()
    f(2) = 0;                       % f-function required for pdepe()
    s(2) = DqDt;                    % s-function required for pdepe()

    % Convert row vectors to column vectors
    c = c'; f = f'; s = s';
end
 
function ic = pde_ic(x)
    ic(1) = 0;                      % Initial concentration at each position, mol/m^3
    ic(2) = 0;                      % Initial adsorbent loading at each position, mol/kg
    ic = ic';
end

function [pl, ql, pr, qr] = pde_bc(xl, ul, xr, ur, t)
    v = 1e-3;                       % Velocity, m/s
    Cin = 100;                      % Concentration at inlet, mol/m^3
 
    % Boundary conditions for the Advection-Diffusion equation
    pl(1) = ul(1) - Cin;            % p-function for pdepe() boundary condition on the left (inlet)
    ql(1) = 0;                      % q-function for pdepe() boundary condition on the left (inlet)
    pr(1) = v*ur(1);                % p-function for pdepe() boundary condition on the right (outlet)
    qr(1) = 1;                      % q-function for pdepe() boundary condition on the right (outlet)

    % Boundary conditions for the Linear Driving Force equation
    pl(2) = 0;                      % p-function for pdepe() boundary condition on the left (inlet)
    ql(2) = 1;                      % q-function for pdepe() boundary condition on the left (inlet)
    pr(2) = 0;                      % p-function for pdepe() boundary condition on the right (outlet)
    qr(2) = 1;                      % q-function for pdepe() boundary condition on the right (outlet)

    % Convert row vectors to column vectors
    pl = pl'; ql = ql'; pr = pr'; qr = qr';
end

% isotherm(c) returns q_star, the adsorbed quantity of fluid under
% concenctration c. This implementation uses the Langmuir adsorption
% isotherm model, which follows the following form:
%
%            q_max * K_L * c
%   q_star = ---------------
%             1 + (K_L * c)
%
% The parameters q_max, the maximum sorbent loading (in this case in units
% of mol/kg) and K_L, the Langmuir coefficient (in this case in units of
% m^3/mol). These two paramaters must be provided and are typically
% determined experimentally. There are other sets of compatible units for
% the components in the above Langmuir model in the concentration basis.
%
% However, if you are provided with model parameters in the partial
% pressure basis, then you will need to convert K_L from a
% pressure-compatible form to a concentration-compatible form. Remember
% that from the ideal gas law, p = c*R*T, where p is the partial pressure
% of a species, c is the concentration of the specieis, R is the specific
% gas constant of the species, and T is the ambient temperature.
%
% The Langmuir adsorption isotherm model is one of many usable adosorption
% isotherm models, and is used in this example for simplicity.
% Additionally, if raw adsorption isotherm data is available, isotherm()
% may be written to directly interpolate between the raw adsorption
% isotherm values. Add a check to ensure that q_star values are negative,
% or else the PDE solution may become unstable.

function q_star = isotherm(c)
    q_max = 1.5e-5;                 % Maximum sorbent loading, mol/kg
    K_L = 1e3;                      % Langmuir coefficient, m^3/mol

    % Evaluate the Langmuir adsorption isotherm model
    q_star = (q_max * K_L * c)/(1 + K_L * c);
end