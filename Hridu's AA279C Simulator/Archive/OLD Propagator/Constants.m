% Constants
clc
clear all
close all

% Heliocentric Reference Frame
state_sun = [0;0;0;0;0;0];

% Asteroid Constants
GM_ast = .0004463; % Placeholder value
R_ast = 8.42; % Placeholder value
inc_ast = 10.82944804139757;

% Satellite Constants
GM_sat = 6.6767e-20 * 2000; % Placeholder value
R_sat = 0.002; % Placeholder value

% Spherical Gravity
order = 34;
degree = 34;
loadgravcoeff
C = zeros(34);
S = zeros(34);
for i = 1:length(GravityCoefficients)
    ord = GravityCoefficients(i, 1);
    deg = GravityCoefficients(i, 2);
    C(ord+1, deg+1) = GravityCoefficients(i, 3);
    S(ord+1, deg+1) = GravityCoefficients(i, 4);
end
clear GravityCoefficients

% Solar Radiation Pressure
Psrp = 4.57e-6;
Csrp_ast = 1; % Placeholder value
Csrp_sat = 1;
Asrp_ast = R_ast^2;
Asrp_sat = R_sat^2;
Msrp_ast = 7.2e15;
Msrp_sat = 2000;

% Drag
A_drag = 0; % Placeholder value
M_drag = 1; % Placeholder value
Cd_drag = 1; % Placeholder value
p0_drag = 1.225e9;
h0_drag = 6378.137;
H_drag = 10;
omega_e = 2*pi/86164;

% Physical Constants
MJD_J2000 = 51544.5;
R_Earth = 6378.137;
R_Sun   = 696000;
R_Moon  = 1738;
AU = 149597870.700;

% Gravitational coefficients
GM_Earth   = 3.986004418e5;                % [km^3/s^2]; WGS-84
GM_Sun     = 1.327124400179870e11;          % [km^3/s^2]; DE405
GM_Moon    = GM_Earth/81.3005600000000044; % [km^3/s^2]; DE405
GM_Mercury = 2.203208048641792e4;          % [km^3/s^2]; DE405
GM_Venus   = 3.248585988264596e5;          % [km^3/s^2]; DE405
GM_Mars    = 4.282831425806710e4;          % [km^3/s^2]; DE405
GM_Jupiter = 1.267127678577960e8;          % [km^3/s^2]; DE405
GM_Saturn  = 3.794062606113726e7;          % [km^3/s^2]; DE405
GM_Uranus  = 5.794549007071872e6;          % [km^3/s^2]; DE405
GM_Neptune = 6.836534063879259e6;          % [km^3/s^2]; DE405
GM_Pluto   = 9.816008877070042e2;          % [km^3/s^2]; DE405
GM_Planets = [GM_Mercury, GM_Venus, GM_Earth,...
    GM_Mars, GM_Jupiter, GM_Saturn, GM_Uranus,...
    GM_Neptune, GM_Pluto];