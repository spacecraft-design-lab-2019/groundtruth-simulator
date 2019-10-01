%EKF test for spin stabalized s/c

%THIS IS NOT GOOD. NEED TO REDO WITH QUATERNION TO AVOID SINGULARITIES

%state: [roll pitch yaw wx wy wz]
%measurements: sun sensor, gyro

dt = 0.05;
tf = 4*60;
t = 0:dt:tf;
n = length(t);
mu = zeros(6,n);
mz = zeros(6,n);
gyro_noise = (pi/180)*(0.1)/(24*60*60);  %accurate to within 0.1 deg/hr
sun_noise = 0.002;
I = [4.056, -0.030, -0.020;    %kg*m2
    -0.030,  4.074, -0.015;
    -0.020, -0.015,  1.897];
sun_eci = [10;0;1]; %random sun direction

%states 
e_i = [0;0;0];
w_i = [0;1;0]*(360/60)*(pi/180);
s_i = [e_i; w_i];
s   = s_i; %assume actual state and estimate state start the same
se  = s_i;
Ee = 0.1*eye(length(s));
Q = 0.01*eye(length(s));
R = sun_noise*eye(6);
R(4,4) = gyro_noise; R(5,5) = gyro_noise; R(6,6) = gyro_noise;

for i =1:n
    %propate actual dynamics (true state)
    a = 10^-9 * sin(t(i)/200);
    torque = [a;0;0];       %some unmodeled torque
    [e_dot,~] = calc_e_dot(s(1:3), s(4:6));
    w_dot = calc_w_dot(s(4:6), torque, I);
    e = s(1:3) + e_dot * dt;
    if(abs(e(3)) > 2*pi)
        e(3) = e(3) - sign(e(3))*2*pi;
    end
    w = s(4:6) + w_dot * dt;
    s = [e;w];
    
    %prediction step (estimate state)
    torque = [0;0;0];       %no torque
    [e_dot,A] = calc_e_dot(se(1:3), se(4:6));
    w_dot = calc_w_dot(se(4:6), torque, I);
    e = se(1:3) + e_dot * dt;
    if(abs(e(3)) > 2*pi)
        e(3) = e(3) - sign(e(3))*2*pi;
    end
    w = se(4:6) + w_dot * dt;
    se = [e;w];
    A2 = eye(6); A2(1:3,1:3) = A;
    Ee = A2*Ee*A2' + Q;
    
    %actual measurements
    gyro_meas = s(4:6) + gyro_noise * randn(3,1);               %GYRO measurement
    [Ry,Rp,Rr] = gen_R_ypr(s(1:3));    
    sun_meas = (Rr*Rp*Ry)'*sun_eci + sun_noise * randn(3,1);    %SUN  measurement
    sun_meas = sun_meas/norm(sun_meas); %normalize, must be direction
    y = zeros(6,1);
    y(1:3) = sun_meas;
    y(4:6) = gyro_meas;
    
    %predicted measurements
    gyro_pred = se(4:6);                %GYRO measurement
    [Ry,Rp,Rr] = gen_R_ypr(se(1:3));    
    sun_pred = (Rr*Rp*Ry)'*sun_eci;     %SUN  measurement
    sun_pred = sun_pred/norm(sun_pred); %normalize, must be direction
    y_hat = zeros(6,1);
    y_hat(1:3) = sun_pred;
    y_hat(4:6) = gyro_pred;
    
    C = gen_C(se, sun_eci);    %measurement matrix coeffecients
    K = Ee*C'*(C*Ee*C'+R)^-1; %Kalman Gain
    se = se+K*(y-y_hat);
    Ee = (eye(6)-K*C)*Ee;
    
    mu(:,i) = se;  %record estimate
    mz(:,i) = s;   %record actual state
end
s_error = mz - mu;

figure, grid on, hold on
plot(t,(180/pi)*s_error(1:3,:))
xlabel('Time [sec]'), ylabel('[deg]')
title('Euler Angle Error')

figure, grid on, hold on
plot(t,(180/pi)*s_error(4:6,:))
xlabel('Time [sec]'), ylabel('[deg/s]')
title('Angular Velocity Error')

figure, grid on, hold on
plot(t,(180/pi)*mu(1:3,:))
plot(t,(180/pi)*mz(1:3,:))
xlabel('Time [sec]'), ylabel('[deg]')
title('Euler Angle Estimate vs Actual')
legend('Re','Pe','Ye','Ra','Pa','Ya');

figure, grid on, hold on
plot(t,(180/pi)*mu(4:6,:))
plot(t,(180/pi)*mz(4:6,:))
xlabel('Time [sec]'), ylabel('[deg/s]')
title('Angular Velocity Estimate vs Actual')
legend('Xe','Ye','Ze','Xa','Ya','Za');

function [e_dot,A] = calc_e_dot( e, w )
%INPUT: euler angles, angular velocity
%OUTPUT: rate of change of euler angles
%propagate actual dynamics
    tan_the = tan(e(2));
    if(abs(tan_the) > 300) %make sure tan(theta) is well defined
        tan_the = 300*sign(tan_the);
        disp('tan x')
    end
    cos_the = cos(e(2));
    if(abs(cos_the) < 10^-4) %make sure cos(theta) is well defined
        cos_the = 10^-4*sign(cos_the);
        disp('cos x')
    end
    A = [1, tan_the*sin(e(1)), tan_the*cos(e(1));
         0, cos(e(1)),              -1*sin(e(1));
         0, sin(e(1))/cos_the,  cos(e(1))/cos_the];
    e_dot = A*w;
end

function w_dot = calc_w_dot(w, tor, I)

w_dot = -I\(Skew(w)*I*w-tor);

% w_dot = zeros(size(w));
% wx = w(1); wy = w(2); wz = w(3);
% Ix = I(1,1); Iy = I(2,2); Iz = I(3,3);
% w_dot(1) = (tor(1) - (Iz-Iy)*wx*wz)/Ix;  %update this to work in body axis
% w_dot(2) = (tor(2) - (Ix-Iz)*wz*wx)/Iy;
% w_dot(3) = (tor(3) - (Iy-Ix)*wx*wy)/Iz;
end

function S = Skew(v)
S = [0 -v(3) v(2)
    v(3) 0 -v(1)
    -v(2) v(1) 0];
end

function [Ry,Rp,Rr] = gen_R_ypr(x)
%given x, [phi, theta, psi]
Ry = [cos(x(3)), -1*sin(x(3)), 0;
      sin(x(3)),    cos(x(3)), 0;
      0,              0,           1];
Rp = [cos(x(2)), 0, sin(x(2));
      0,           1,     0;
   -1*sin(x(2)), 0, cos(x(2))];
Rr = [1,0,0;
      0, cos(x(1)), -1*sin(x(1));
      0, sin(x(1)),    cos(x(1))];
end

function C = gen_C(se, sun_eci)
%given state estimate, expected sun position
%return measurement sensitivity matrix
%measurements are sun body xyz; gyro body xyz
C = eye(6);
ca = cos(se(1)); sa = sin(se(1));
cb = cos(se(2)); sb = sin(se(2));
cc = cos(se(3)); sc = sin(se(3));
sx = sun_eci(1); sy = sun_eci(2); 
sz = sun_eci(3);
xa = (0)*sx+(ca*sb*cc-sa*sc)*sy+(ca*sc+sa*sb*cc)*sz;
xb = (-1*sb*cc)*sx+(sa*cb*cc)*sy+(-1*ca*cb*cc)*sz;
xc = (-1*cb*sc)*sx+(-1*sa*sb*sc+ca*cc)*sy+(sa*cc+ca*sb*sc)*sz;
ya = (0)*sx+(-1*sa*cc-ca*sb*sc)*sy+(ca*cc-sa*sb*sc)*sz;
yb = (sb*sc)*sx+(-1*sa*cb*sc)*sy+(ca*cb*sc)*sz;
yc = (-1*cb*cc)*sx+(-1*ca*sc-sa*sb*cc)*sy+(-1*sa*sc+ca*sb*cc)*sz;
za = (0)*sx+(-1*ca*cb)*sy+(-1*sa*cb)*sz;
zb = (cb)*sx+(sa*sb)*sy+(-1*ca*sb)*sz;
zc = 0;

C(1:3, 1:3) = zeros(3,3);
%C(1:3, 1:3) = [xa, xb, xc;
%               ya, yb, yc;
%               za, zb, zc]; 
end
