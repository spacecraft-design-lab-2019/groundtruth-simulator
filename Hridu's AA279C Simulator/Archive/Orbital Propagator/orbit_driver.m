%--------------------------------------------------------------------------
% SCRIPT: orbit_driver.m
% AUTHOR: Adam Zufall & Hridu Jain
% DATE:   April 23, 2019

% NOTES: (none)
%--------------------------------------------------------------------------

% clc; clear all; close all;

%% Initialize Orbit

GM_Moon = 4.9048695e3;      % [km^3/s^2]
r_Moon = 1737.4;            % [km] - Volumetric Mean Radius

r_p = r_Moon + 50;          % [km]
r_a = r_Moon + 125;         % [km]

a = .5*(r_p+r_a);           % [km]
e = (r_a-r_p)/(r_a+r_p);
i = 0;                      % [rad]
Omega = 0;                  % [rad]
omega = 0;                  % [rad] % NOTE UPDATE THIS VALUE
nu = 0;                     % [rad]

state_initial = calcPositionAndVelocity(a, e, i, Omega, omega, nu, GM_Moon);


%% Propogate
t_span = [0, 5*60]; % sim time = 5 min
options = odeset('RelTol',1e-6);
[t_out, state_out] = ode45(@(t, state) f(t, state, GM_Moon), t_span, state_initial, options);


%% Plot

% % 2-D orbit around Moon-Centered Inertial Coordinate System
% figure, grid on, hold on
% plot(state_out(:, 1), state_out(:, 2))



%% State Dot
function statedot = f(t, state, GM)
    statedot = zeros(6,1);
    statedot(1:3) = state(4:6);
    statedot(4:6) = -GM*state(1:3)./norm(state(1:3))^3;
end