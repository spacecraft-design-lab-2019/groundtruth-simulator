%EKF test for spin stabalized s/c

%REDO WITH QUATERNION TO AVOID SINGULARITIES

%state: [q1 q2 q3 q4 wx wy wz]
%measurements: sun sensor, gyro

dt = 0.05;
tf = 10*60*60;
t = 0:dt:tf;
n = length(t);
mu = zeros(7,n);
mz = zeros(7,n);
gyro_noise = (pi/180)*(0.1)/(24*60*60);  %accurate to within 0.1 deg/hr
sun_noise = 0.002;
I = [4.056, -0.030, -0.020;    %kg*m2
    -0.030,  4.074, -0.015;
    -0.020, -0.015,  1.897];
sun_eci = [10;0;1]; %random sun direction

%states 
q_i = [-0.145392350642603;0;0;0.989374077068233];
w_i = [0;0;1]*(360/60)*(pi/180);
s_i = [q_i; w_i];
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
    [q_dot,~] = calc_q_dot(s(1:4), s(5:7));
    w_dot = calc_w_dot(s(5:7), torque, I);
    q = s(1:4) + q_dot * dt;
    q = q/norm(q);
    w = s(5:7) + w_dot * dt;
    s = [q;w];
    
    %prediction step (estimate state)
    torque = [0;0;0];       %no torque
    [q_dot,A] = calc_q_dot(se(1:4), se(5:7));
    w_dot = calc_w_dot(se(5:7), torque, I);
    q = se(1:4) + q_dot * dt;
    q = q/norm(q);
    w = se(5:7) + w_dot * dt;
    se = [q;w];
    A2 = eye(7); A2(1:4,1:4) = A;
    Ee = A2*Ee*A2' + Q;
    
    %actual measurements
    gyro_meas = s(5:7) + gyro_noise * randn(3,1);               %GYRO measurement
    Rot = calc_R(s(1:4));
    sun_meas = Rot*sun_eci + sun_noise * randn(3,1);    %SUN  measurement
    sun_meas = sun_meas/norm(sun_meas); %normalize, must be direction
    y = zeros(6,1);
    y(1:3) = sun_meas;
    y(4:6) = gyro_meas;
    
    %predicted measurements
    gyro_pred = se(5:7);                %GYRO measurement
    Rot = calc_R(se(1:4));    
    sun_pred = Rot*sun_eci;     %SUN  measurement
    sun_pred = sun_pred/norm(sun_pred); %normalize, must be direction
    y_hat = zeros(6,1);
    y_hat(1:3) = sun_pred;
    y_hat(4:6) = gyro_pred;
    
    C = gen_C(se, sun_eci);    %measurement matrix coeffecients
    K = Ee*C'*(C*Ee*C'+R)^-1; %Kalman Gain
    se = se+K*(y-y_hat);
    Ee = (eye(7)-K*C)*Ee;
    
    mu(:,i) = se;  %record estimate
    mz(:,i) = s;   %record actual state
end
s_error = mz - mu;

figure, grid on, hold on
plot(t,s_error(1:4,:))
xlabel('Time [sec]'), ylabel('[deg]')
title('Quaternion Error')

figure, grid on, hold on
plot(t,(180/pi)*s_error(5:7,:))
xlabel('Time [sec]'), ylabel('[deg/s]')
title('Angular Velocity Error')

figure, grid on, hold on
plot(t,mu(1:4,:))
plot(t,mz(1:4,:))
xlabel('Time [sec]'), ylabel('[deg]')
title('Euler Angle Estimate vs Actual')
legend('q1e','q2e','q3e','q4e','q1a','q2a','q3a','q4a');

figure, grid on, hold on
plot(t,(180/pi)*mu(5:7,:))
plot(t,(180/pi)*mz(5:7,:))
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
sun_eci = sun_eci / norm(sun_eci);
C = zeros(6,7);
q1 = se(1); q2 = se(2);
q3 = se(3); q4 = se(4);

temp = [ q1,  q2,  q3;
        -q2,  q1, -q4;
        -q3,  q4,  q1;
         q4,  q3, -q3;
         q2, -q1,  q4;
         q1,  q2,  q3;
        -q4, -q3,  q2;
        -q3,  q4,  q1;
         q3, -q4, -q1;
         q4,  q3, -q2;
         q1,  q2,  q3;
         q2, -q1,  q4];
temp = 2*temp*sun_eci;
x1 = temp(1); x2 = temp( 2); x3 = temp( 3); x4 = temp( 4);
y1 = temp(5); y2 = temp( 6); y3 = temp( 7); y4 = temp( 8);
z1 = temp(9); z2 = temp(10); z3 = temp(11); z4 = temp(12);
C(1:3, 1:4) = [ x1, x2, x3, x4;
                y1, y2, y3, y4;
                z1, z2, z3, z4];
C(1:3, 1:4) = zeros(3,4);
C(4:6, 5:7) = eye(3);
end

function [q_dot, A] = calc_q_dot(q, w)

wx = w(1); wy = w(2); wz = w(3);
OMEGA = [   0,    wz, -1*wy, wx;
        -1*wz,     0,    wx, wy;
           wy, -1*wx,     0, wz;
        -1*wx, -1*wy, -1*wz,  0];
q_dot = 0.5*OMEGA*q;
A = 0.5*OMEGA;
end

function R_eci2body = calc_R(q)
eta=q(4);
e1=q(1);
e2=q(2);
e3=q(3);

R_eci2body = [e1^2-e2^2-e3^2+eta^2, 2*(e1*e2+eta*e3),        2*(e1*e3-eta*e2);
              2*(e1*e2-eta*e3),     -1*e1^2+e2^2-e3^2+eta^2, 2*(e2*e3+eta*e1);
              2*(e1*e3+eta*e2),     2*(e2*e3-eta*e1),        -1*e1^2-e2^2+e3^2+eta^2];
          
end