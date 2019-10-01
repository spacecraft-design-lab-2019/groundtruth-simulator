%--------------------------------------------------------------------------
% SCRIPT: driverADCS.m
% AUTHOR: Adam Zufall & Hridu Jain
% DATE:   April 30, 2019

% NOTES:
% driver for ADCS.slx
%--------------------------------------------------------------------------

clc; clear all; close all;

%% Global Constants
GM_Moon = 4.9048695e3;      % [km^3/s^2]
r_Moon = 1737.4;            % [km] - Volumetric Mean Radius


%% Initialize Orbit

r_p = r_Moon + 50;          % [km]
r_a = r_Moon + 125;         % [km]

a = .5*(r_p+r_a);           % [km]
e = (r_a-r_p)/(r_a+r_p);
i = 0;                      % [rad]
Omega = 0;                  % [rad]
omega = 0;                  % [rad] % NOTE UPDATE THIS VALUE
nu = 0;                     % [rad]

oe_initial = [a, e, i, Omega, omega, nu];


%% Run

sim('ADCS.slx')