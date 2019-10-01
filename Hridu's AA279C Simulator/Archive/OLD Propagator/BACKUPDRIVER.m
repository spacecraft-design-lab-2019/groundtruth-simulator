% ----------------------------------------------------
% BACKUP DRIVER


%SLAB 2017 REU
% Written by: Hridu Jain
% Date: August 24, 2017
% ----------------------------------------------------

clc
clear all
% close all

X  = load('constants.mat');

jd_start = datenum(2017, 08, 01, 00, 00, 00) + 1721058.5;
jd_end   = datenum(2018, 02, 01, 00, 00, 00) + 1721058.5;
stepsize = .5;

% Initial Parameters of Asteroid w.r.t. Sun
xin = -1.831318104214490E+08;
yin = -6.194513053903891E+07;
zin = -3.560695348353375E+07;
vxin = 3.684009600658319E+00;
vyin = -2.688513550415339E+01;
vzin = -2.317567534994446E+00;
state_ast_initial = [xin; yin; zin; vxin; vyin; vzin];

% Initial Parameters of Satellite w.r.t. Asteroid
a = 1000;               % [km]
e = 0;
i = 0;                  % [rad]
Omega = 0;              % [rad]
omega = 0;              % [rad]
E = 0;                  % [rad]

% [r_initial, v_initial] = calcPositionAndVelocity(a, e, i, Omega, omega, E, X.GM_ast);
% state_sat_initial = ECItoSCI([r_initial + state_ast_initial(1:3); v_initial], X.inc_ast);
state_sat_initial = state_ast_initial;
state_sat_initial(1) = state_sat_initial(1) + 1000;

sim('Propagator.slx')

% Calculating Relative Orbital Elements
roe = zeros(length(oe_sat), 6);
for i = 1:length(oe_sat)
    roe(i, :) = oe_ast(i, 1)*oe2roe(oe_sat(i, :), oe_ast(i, :));
end
tout = tout - tout(1);
plotroe(tout./643, roe);

