load PS4_info_1st20.mat
load initial_trans.mat

t = tnew;
tstep = t(2)-t(1);

R_earth = NormREarth;
R_sun = NormRSun;
mu_sun = 1.32712440018e11;
mu_earth = 3.986004418e5;
Ix = 184.9343576; % kg * m^2
Iy = 321.8711548; % kg * m^2
Iz = 423.2775247; % kg * m^2
I = [Ix 0 0;
    0 Iy 0;
    0 0 Iz];
C = initial_trans;
beta4sq = 1/4 * (1 + trace(C));
beta1sq = 1/4 * (1 + 2*C(1,1) - trace(C));
beta2sq = 1/4 * (1 + 2*C(2,2) - trace(C));
beta3sq = 1/4 * (1 + 2*C(3,3) - trace(C));
beta1 = sqrt(beta1sq);
beta2 = (C(3,1) + C(1,3))/(4*beta1);
beta3 = (C(1,2) + C(2,1))/(4*beta1);
beta4 = (C(2,3) - C(3,2))/(4*beta1);
% q0 = [0; 0; 0; 1];
q0 = [beta1; beta2; beta3; beta4];
wx0 = .5;
wy0 = .2;
wz0 = .4;
w0 = [wx0; wy0; wz0]; % rad/s

state0 = [Ix; Iy; Iz; w0; q0];


n = length(t);

c_sun = [cxsun, cysun, czsun];
M_sun = GravTorques(mu_sun,R_sun,c_sun,I);

c_earth = [cxearth, cyearth, czearth];
M_earth = GravTorques(mu_earth,R_earth,c_earth,I);