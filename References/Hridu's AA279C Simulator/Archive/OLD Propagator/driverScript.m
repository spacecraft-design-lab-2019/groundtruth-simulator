% ----------------------------------------------------
% SLAB 2017 REU
% Written by: Hridu Jain
% Date: August 24, 2017
% ----------------------------------------------------

clc
clear all
% close all

X  = load('constants.mat');

jd_start = datenum(2017, 08, 01, 00, 00, 00) + 1721058.5;
%jd_end   = datenum(2017, 08, 05, 00, 00, 00) + 1721058.5;
jd_end = jd_start + 0.4*365.25;
stepsize = (jd_end - jd_start)/600;

jd_start = jd_start*24*60*60;
jd_end   = jd_end *24*60*60;
stepsize = (jd_end - jd_start)/2000;


% Initial Parameters of Asteroid w.r.t. Sun
xin = -1.831318104214490E+08;
yin = -6.194513053903891E+07;
zin = -3.560695348353375E+07;
vxin = 3.684009600658319E+00;
vyin = -2.688513550415339E+01;
vzin = -2.317567534994446E+00;
state_ast_initial = [xin; yin; zin; vxin; vyin; vzin];

% Initial Parameters of Satellite w.r.t. Asteroid
a = 100;               % [km]
e = 0;
i = 0;                  % [rad]
Omega = 0;              % [rad]
omega = 0;              % [rad]
E = 0;                  % [rad]

% [r_initial, v_initial] = calcPositionAndVelocity(a, e, i, Omega, omega, E, X.GM_ast);
% state_sat_initial = ECItoSCI([r_initial + state_ast_initial(1:3); v_initial], X.inc_ast);
xin2 = -1.828016060904607E+08;
yin2 = -6.426392646769845E+07;
zin2 = -3.580487053967395E+07;
vxin2 = 3.959196098092854E+00;
vyin2 = -2.679022882188029E+01;
vzin2 = -2.263860639369534E+00;
state_sat_initial = [xin2; yin2; zin2; vxin2; vyin2; vzin2];

xin = -1.841318104214490E+08; %let's try something different
yin = -6.194513053903891E+07;
zin = -3.560695348353375E+07;
vxin = 3.684009600658319E+00;
vyin = -2.708513550415339E+01;
vzin = -2.317567534994446E+00;
state_sat_initial = [xin2; yin2; zin2; vxin2; vyin2; vzin2];

sim('Propagator.slx')

% Calculating Relative Orbital Elements
roe = zeros(length(oe_sat), 6);
for i = 1:length(oe_sat)
    roe(i, :) = oe_ast(i, 1)*oe2roe(oe_sat(i, :), oe_ast(i, :));
end
tout = tout - tout(1);
plotroe(tout./643, roe);


%% plot s/c and asteroid's relative orbits
r_p1 = squeeze(state_planets(1,1:3,:));
r_p2 = squeeze(state_planets(2,1:3,:));
r_p3 = squeeze(state_planets(3,1:3,:));
r_p4 = squeeze(state_planets(4,1:3,:));

figure
hold on
plot3([0],[0],[0],'o')
plot3(state_sat(:,1), state_sat(:,2), state_sat(:,3))
plot3(state_ast(:,1), state_ast(:,2), state_ast(:,3))
plot3(r_p1(1,:), r_p1(2,:), r_p1(3,:))
plot3(r_p2(1,:), r_p2(2,:), r_p2(3,:))
plot3(r_p3(1,:), r_p3(2,:), r_p3(3,:))
plot3(r_p4(1,:), r_p4(2,:), r_p4(3,:))
legend('Sun','S/C','Eros','Mercury','Venus','Earth','Mars')
axis('equal')

figure
plot3(state_sat(:,1)-state_ast(:,1), state_sat(:,2)-state_ast(:,2), state_sat(:,3)-state_ast(:,3))
title('S/C Relative to Eros')