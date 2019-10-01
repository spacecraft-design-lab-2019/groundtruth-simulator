function [Kp, Kd] = calc_gains(I,zeta, rise_time)
%given inertia matrix, damping ratio, rise time (seconds)
%return Kp and Kd for system
Ixx = I(1,1);
Iyy = I(2,2);
Izz = I(3,3);

gain_margin = 1.5;

w_n = (pi-atan(sqrt(1-zeta^2)/zeta))/(rise_time*sqrt(1-zeta^2)); %natural frequency

Kpx = gain_margin * Ixx * w_n^2;
Kpy = gain_margin * Iyy * w_n^2;
Kpz = gain_margin * Izz * w_n^2;
Kp = [Kpx; Kpy; Kpz];

Kdx = 2*gain_margin*zeta*Ixx*w_n;
Kdy = 2*gain_margin*zeta*Iyy*w_n;
Kdz = 2*gain_margin*zeta*Izz*w_n;
Kd = [Kdx; Kdy; Kdz];

end

