%init script for sim
clc, clear all, close all

%% simulation parameters
delta_t  = 0.2;   %seconds
end_time = 5*60;  %seconds

%% initial state
SC.initial_state.r_eci = [6380000; 0; 0];       %m
SC.initial_state.v_eci = [0; 8000; 0];          %m/s
SC.initial_state.q     = [1; 0; 0; 0];
SC.initial_state.e     = quat2eul(SC.initial_state.q);
SC.initial_state.w     = (pi/180) * [0; 0; 15]; %rad/s
SC.initial_state.date  = [2019, 05, 06, 12, 00, 00];

%% mass properties
SC.mass_properties.inertia = [4.056, -0.030, -0.020;    %kg*m2
                             -0.030,  4.074, -0.015;
                             -0.020, -0.015,  1.897];
SC.mass_properties.mass = 54;       %kg, dry mass for now
SC.mass_properties.cg_offset = [0.33, 0.24, 0.16]/100;  %m
SC.mass_properties.R_body2principal = [1,0,0;
                                       0,1,0;
                                       0,0,1];

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
%% RCS
SC.RCS.dir = [1, 0, 0,   1, 0, 0,  -1, 0, 0,  -1, 0, 0;
              0, 1, 0,   0,-1, 0,   0, 1, 0,   0,-1, 0;
              0, 0, 1,   0, 0,-1,   0, 0,-1,   0, 0, 1];
SC.RCS.moment_arms = [0, 3, 3,   0, 3, 3,   0,-3,-3,   0,-3,-3;        %m
                      3, 0, 3,  -3, 0,-3,   3, 0, 3,  -3, 0,-3;
                      4, 4, 0,  -4,-4, 0,  -4,-4, 0,   4, 4, 0]/10;
SC.RCS.max_thrust = 0.1;  %N
%% sensors
SC.sensors.sun.noise = 0.002;                       %accurate to 0.2 degrees
SC.sensors.gyro.bias = [0;0;0];                     %rad/s 
SC.sensors.gyro.noise = (pi/180)*(0.1)/(24*60*60);  %accurate to within 0.1 deg/hr
SC.sensors.star.noise = 20/(60*60);                 %accurate to 20 arc seconds
SC.sensors.star.w_limit = (pi/180)*2;               %star tracker cannot operate at over 2 deg/sec

%% constants
SS.constants.GM_Earth = 3.9860044188   * 10^14; %m3/s2
SS.constants.GM_Moon  = 4.90486959     * 10^12; %m3/s2
SS.constants.GM_Sun   = 1.327124400189 * 10^20; %m3/s2
SS.constants.igrf = load('igrf_constants.mat');
%SS.constants.igrf.yearLastEpoch;
%SS.constants.igrf.lastgh;
%SS.constants.igrf.ghslope;