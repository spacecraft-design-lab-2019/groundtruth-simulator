%init script for sim
clc, clear all, close all

%% simulation parameters
delta_t  = 0.1;   %seconds
end_time = 30*60;  %seconds

%% initial state
%SC.initial_state.r_eci = [-6293110.28347100;-2087790.56184149;-1457505.19446010];   %ISS orbit
SC.initial_state.r_eci = [63684080.5957580;348960520.146146;133767644.928460];      %moon orbit
%SC.initial_state.v_eci = [330.999339309726; -5019.55238831980;5783.16246153759];    %ISS orbit
SC.initial_state.v_eci = [-1027.49637842178;1766.36877974272;-371.263903666297];    %moon orbit

%SC.initial_state.r_mci = [1737150; 0;0 ];
%SC.initial_state.v_mci = [0; 1680.4; 0];

SC.initial_state.q     = [-0.145392350642603;0;0;0.989374077068233];
SC.initial_state.e     = [-0.291819050933452; 0; 0];    %quat2eul(SC.initial_state.q);
% SC.initial_state.w     = [0;0;7]*pi/180;                %RPM -> rad/s
SC.initial_state.w     = [1;5;2Usin]*pi/180;
SC.initial_state.date  = [2019, 05, 06, 12, 00, 00];    %JD = 2458610.00000
SC.initial_state.date  = [2019, 05, 07, 16, 59, 39.9];  %JD = 2458611.208101

%% mass properties
SC.mass_properties.inertia = [4.056, -0.030, -0.020;    %kg*m2
                             -0.030,  4.074, -0.015;
                             -0.020, -0.015,  1.897];
SC.mass_properties.inertia = [1.8965, 0, 0;    %kg*m2
                             0,  4.0963, 0;
                             0, 0,  4.0963];
SC.mass_properties.inertia = [0.8965, 0, 0;    %kg*m2
                             0,  2.0963, 0;
                             0, 0,  4.0963];
SC.mass_properties.mass = 54;                           %kg, dry mass for now
SC.mass_properties.cg_offset = [0.33, 0.24, 0.16]/100;  %m
SC.mass_properties.R_body2principal = [1,0,0;
                                       0,1,0;
                                       0,0,1];

%% controller gains and time constants
SC.ctrl.zeta = 0.707;
SC.ctrl.rise_time = 20; %seconds
[Kp, Kd] = calc_gains(SC.mass_properties.inertia,SC.ctrl.zeta, SC.ctrl.rise_time);

SC.ax3.gains = [Kp; Kd];
SC.ax3.time_constant = 0.1;
SC.spin.gains = [0*Kp; Kd];
SC.spin.time_constant = 0.1;
SC.spin_up.gains = [Kp; Kd];
SC.spin_up.time_constant = 0.1;

%% structure
SC.structure.unitSurfaces = [1,-1,0, 0,0, 0;
                             0, 0,1,-1,0, 0;
                             0, 0,0, 0,1,-1];
SC.structure.area_per_face = [0.48, 0.48, 0.48, 0.48, 0.36, 0.36]; %m2 %+x,-x,+y,-y,+z,-z
SC.structure.moment_arms = [0,        0,        0.001629,-0.001623,-0.002443, 0.002443; %moment arms for SRP and drag
                           -0.001623, 0.001629, 0,        0,        0.003257,-0.003257;
                            0.002443,-0.002443,-0.003257, 0.003257, 0,        0];
SC.structure.CD = 1.5;  %coeffecient of drag
SC.structure.residual_dipole = [0;0;0]; %Am2

SC.structure.srp_coeff.absorption   = 0.79;
SC.structure.srp_coeff.spec_reflec  = 0.21;
SC.structure.srp_coeff.diff_reflec  = 0.00;

%% RCS
SC.RCS.dir = [1, 0, 0,   1, 0, 0,  -1, 0, 0,  -1, 0, 0;
              0, 1, 0,   0,-1, 0,   0, 1, 0,   0,-1, 0;
              0, 0, 1,   0, 0,-1,   0, 0,-1,   0, 0, 1];
SC.RCS.moment_arms = [0, 3, 3,   0, 3, 3,   0,-3,-3,   0,-3,-3;        %m
                      3, 0, 3,  -3, 0,-3,   3, 0, 3,  -3, 0,-3;
                      4, 4, 0,  -4,-4, 0,  -4,-4, 0,   4, 4, 0]/10;
SC.RCS.max_thrust = 0.1/10;  %N
SC.RCS.num_bits = 256;
SC.RCS.deadband = 1/25000; %cmd, unitless
SC.RCS.bias = 0*(pi/180) * 0.05 * randn(12,3); %misalignment, rad
SC.RCS.noise = 0*0.01;
%% sensors
%software
SC.sensors.gyro.gains = [1;1;1];
SC.sensors.gyro.offset = [0;0;0];
SC.sensors.gyro.rot_mat = [1,0,0; 0,1,0; 0,0,1];
SC.sensors.gyro.wc = 600 * (pi/180);

SC.sensors.sun.pre_wc = 600 * (pi/180);
SC.sensors.sun.gains = [1;1;1;1;1;1];
SC.sensors.sun.offset = [0;0;0;0;0;0];
SC.sensors.sun.threshold = (1/10);
SC.sensors.sun.dir_wc = 600 * (pi/180);

SC.sensors.star.wc = 20 * (pi/180);

SC.sensors.weights = [1;1;1;0.1];

%hardware
noise = 1;
if(noise)
    SC.sensors.sun.noise = 0.002;                       %accurate to 0.2 degrees
    SC.sensors.gyro.bias = [0;0;0];                     %rad/s 
    SC.sensors.gyro.noise = (pi/180)*(0.1)/(24*60*60);  %accurate to within 0.1 deg/hr
    SC.sensors.star.noise =  5/(60*60);                 %accurate to 20 arc seconds
    SC.sensors.star.w_limit = (pi/180)*2;               %star tracker cannot operate at over 2 deg/sec
    temp_angle = (pi/180) * (20/3600) * randn(1,3); %misalignment, rad
    ax = temp_angle(1); ay = temp_angle(2); az = temp_angle(3);
    R_bias = [    1,    az, -1*ay;
              -1*az,     1,    ax;
                 ay, -1*ax,     1];
    q_bias = dcm2quat(R_bias); %scalar first for MATLAB function  
    SC.sensors.star.bias = q_bias/norm(q_bias);                     %misalignment & bias
else
    SC.sensors.sun.noise = 0.00;                       %accurate to 0.2 degrees
    SC.sensors.gyro.bias = [0;0;0];                     %rad/s 
    SC.sensors.gyro.noise = (pi/180)*(0.0)/(24*60*60);  %accurate to within 0.1 deg/hr
    SC.sensors.star.noise = 0/(60*60);                 %accurate to 20 arc seconds
    SC.sensors.star.w_limit = (pi/180)*2;               %star tracker cannot operate at over 2 deg/sec
    SC.sensors.star.bias = [1;0;0;0];                     %misalignment / bias
    SC.sensors.sun.dir_wc = 1000*600 * (pi/180);
    SC.sensors.star.wc = 1000*20 * (pi/180);
    SC.sensors.gyro.wc = 1000*600 * (pi/180);
end

%% EKF
R = 5*SC.sensors.sun.noise*eye(6);
R(4:6,4:6) = SC.sensors.gyro.noise*eye(3); 
SC.estimator.sensor_noise = R;
R = 5*SC.sensors.sun.noise*eye(10);
R(4:6,4:6) = SC.sensors.gyro.noise*eye(3); 
R(7:10, 7:10) = SC.sensors.star.noise*eye(4);
SC.estimator.sensor_noise_w_star = R;
Q = 10^-7*eye(7);
Q(5:7,5:7) = 10^-4*eye(3); 
SC.estimator.process_noise = Q; 
P = SC.sensors.star.noise * eye(7);
P(5:7, 5:7) = SC.sensors.gyro.noise * eye(3);
SC.estimator.initial_sigma = P;
clear R, clear P, clear Q;
%% constants
SS.constants.GM_Earth = 3.9860044188   * 10^14; %m3/s2
SS.constants.GM_Moon  = 4.90486959     * 10^12; %m3/s2
SS.constants.GM_Sun   = 1.327124400189 * 10^20; %m3/s2
SS.constants.igrf = load('igrf_constants.mat');
%SS.constants.igrf.yearLastEpoch;
%SS.constants.igrf.lastgh;
%SS.constants.igrf.ghslope;
SS.constants.Obliquity_Earth = 23.43929111;     %degrees
SS.constants.Obliquity_Moon  =  6.68;           %degrees
SS.constants.orbit_coef = load('sun_and_moon_Aug2018_Sept2021.mat');
SS.constants.radius_moon  = 1737100;    %m
SS.constants.radius_earth = 6371000;    %m


%% RUN

sim('ADCS_r2017b.slx')
plotting_script

