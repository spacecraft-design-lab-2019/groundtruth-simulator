%--------------------------------------------------------------------------
% SCRIPT: dynamic_first_try.m
% AUTHOR: Adam Zufall
% DATE:   April 9, 2019

% NOTES:
% First Dynamic Test
% help determine mass properties necessary to achieve pointing requirements
%--------------------------------------------------------------------------

clc; clear all; close all;

%% Initialize
e_o = [0; 0; 0];            % initial Euler angles, rad
q_o = [1; 0; 0; 0];         % initial quaternion
w_o = (pi/180)*[5;-3;7];     % initial rotation rate, rad/s
s_o = [q_o; w_o];           % initial state, q w

%% Moment of Inertia
% play around with different moments of inertia
% CG = [3, -5, 2]/100000;       %center of gravity location, m
% I = calc_I(CG);             %moment of inertia, kg*m2

I = 10^0 * [6.57, 0, 0;
            0, 6.57, 0;
            0, 0, 12.22];
Ix = I(1,1); Iy = I(2,2); Iz = I(3,3);
% ant_pt = [0.5, 0.0, 0.0;  %arbitrary points of interest
%           0.0, 0.0, 0.5;
%          -0.5, 0.0, 0.0;
%           0.0, 0.0,-0.5]; %where the antenae meet the body, m

%pt_dir = [0;0;1];

%% Attitude Propogation
dt = 0.01;
t_o = 0;
t_f = 5*60; % sim time is 30 minutes
t = t_o:dt:t_f;
n = length(t);
s = zeros(7,n);
s(:, 1) = s_o;

for i=2:n
    %euler, first order integration
    %s_dot = rot_dyn(s(:,i-1));
    %s(:,i) = s(:,i-1)+s_dot*dt;  %Euler method
    
    %rk4, higher order integration
    k1 = rot_dyn(s(:,i-1), I);
    k2 = rot_dyn(s(:,i-1)+0.5*dt*k1, I);
    k3 = rot_dyn(s(:,i-1)+0.5*dt*k2, I);
    k4 = rot_dyn(s(:,i-1)+    dt*k3, I);
    s(:,i) = s(:,i-1) + dt*(k1+2*k2+2*k3+k4)/6;
    
%     % wrap as needed - for EULER ANGLE representation
%     if(abs(s(1,i)) > pi)
%         s(1,i) = s(1,i)-2*pi*sign(s(1,i));
%     end
%     if(abs(s(3,i)) > pi)
%         s(3,i) = s(3,i)-2*pi*sign(s(3,i));
%     end

    %make sure q stays a unit vector - for QUATERNION representation
    s(1:4, i) = s(1:4, i)/norm(s(1:4, i));

    %apply fictional control, nutation damping
    %s(4:5,i) = 0.9999*s(4:5,i);
end

wx = s(5,:);
wy = s(6,:);
wz = s(7,:);

a_wxy = sqrt(w_o(1)^2 + w_o(2)^2);
phi = atan2(w_o(2), w_o(1));
lam = (Iz-Ix)/Ix; %(2*pi)/w_o(3);
for i=1:n
    a_wx(i) = a_wxy*(cos(lam*t(i)+phi));
    a_wy(i) = a_wxy*(sin(lam*t(i)+phi));
    
    error(i) = norm([a_wx(i)-wx(i); a_wy(i)-wy(i)]);
end

%determine interest points motion
% ant_pts = zeros(4,3,n);
% ant_ang = zeros(4,  n);
% pt_dirs = zeros(3,  n);
% pt_err  = zeros(1,  n);
% for i=2:n
%     R = calc_R(s(1:3,i));
%     ant_pts(1,:,i) = R*ant_pt(1,:)';
%     ant_pts(2,:,i) = R*ant_pt(2,:)';
%     ant_pts(3,:,i) = R*ant_pt(3,:)';
%     ant_pts(4,:,i) = R*ant_pt(4,:)';
%     pt_dirs(  :,i) = R*pt_dir / norm(R*pt_dir);
%     
%     ant_ang(1,i) = asin(ant_pts(1,3,i));
%     ant_ang(2,i) = asin(ant_pts(2,3,i));
%     ant_ang(3,i) = asin(ant_pts(3,3,i));
%     ant_ang(4,i) = asin(ant_pts(4,3,i));
%     pt_err(   i) = acos(dot(pt_dirs(:,i), pt_dir));
% end

%% calculate momentum and energy ellipsoids

for i=1:n
    twoT(i) = wx(i)*wx(i)*I(1,1)        + wy(i)*wy(i)*I(2,2)        + wz(i)*wz(i)*I(3,3); 
    Lsq(i)  = wx(i)*wx(i)*I(1,1)*I(1,1) + wy(i)*wy(i)*I(2,2)*I(2,2) + wz(i)*wz(i)*I(3,3)*I(3,3); 
end
num = 50;
rxe = sqrt(mean(twoT) / I(1,1));
rye = sqrt(mean(twoT) / I(2,2));
rze = sqrt(mean(twoT) / I(3,3));
[xe, ye, ze] = ellipsoid(0,0,0, rxe, rye, rze, num);

L = sqrt(mean(Lsq));
rxm = L / I(1,1);
rym = L / I(2,2);
rzm = L / I(3,3);
[xm, ym, zm] = ellipsoid(0,0,0, rxm, rym, rzm, num);

rxp = ((I(2,2)-Lsq(1)*twoT(1))*(I(3,3)-Lsq(1)*twoT(1)))^-1;
ryp = ((I(1,1)-Lsq(1)*twoT(1))*(I(3,3)-Lsq(1)*twoT(1)))^-1;
rzp = ((I(2,2)-Lsq(1)*twoT(1))*(I(1,1)-Lsq(1)*twoT(1)))^-1;
[xp, yp, zp] = ellipsoid(0,0,0, rxp, ryp, rzp, num);

omx = linspace(min(wx), max(wx), num);
omy = linspace(min(wy), max(wy), num);

L2 = mean(Lsq);
T2 = mean(twoT);
for a = 1:num
    for b = 1:num
        p1 = -1*(Ix-L2/T2)*Ix*omx(a)*omx(a);
        p2 = -1*(Iy-L2/T2)*Iy*omy(b)*omy(b);
        p3 = (Iz-L2/T2)*Iz;
        omz(a,b) = sqrt((p1+p2)/p3);
    end
end

%% calculate euler angles from quaternion
for i=1:n
    [rol(i), pit(i), yaw(i)] = quat2eul(s(1:4,i));
end
eul_ang = [rol; pit; yaw];

%% plot
figure, grid on, hold on
%plot(t, s(4:6,:)*180/pi)
plot(t, s(5:7,:)*180/pi, 'LineWidth',2)
xlabel('Time, seconds')
ylabel('\omega, deg/sec')
title('Rotation Rate')
legend('X','Y','Z')

figure, grid on, hold on
%plot3(s(4,:), s(5,:), s(6,:))
plot3(s(5,:), s(6,:), s(7,:))
title('Angular Velocity Profile')
xlabel('X'), ylabel('Y'), zlabel('Z')

figure, grid on, hold on
plot(t, eul_ang(1:3,:)*180/pi)
title('Euler Angles')
xlabel('time, seconds'), ylabel('angle, degrees')
legend('Roll','Pitch','Yaw')

figure, grid on, hold on
plot(t, s(1:4,:)*180/pi)
title('Quaternion Components')
xlabel('time, seconds'), ylabel('angle, degrees')
legend('q1','q2','q3','q4')

% figure, grid on, hold on
% plot(t, 1000*squeeze(ant_pts(1,3,:)))
% title('Antenae 1')
% xlabel('Time, seconds'), ylabel('millimeters')
% %legend('X','Y','Z')
% 
% figure, grid on, hold on
% plot(t, ant_ang(1,:)*180/pi)
% xlabel('Time'), ylabel('Angle, degrees')
% title('Antenae 1')

% figure, grid on, hold on
% plot(t, pt_err*180/pi)
% xlabel('Time, seconds'), ylabel('Angle, degrees')
% title('Pointing Error')

figure, grid on, hold on
plot(t, twoT)
xlabel('Time, seconds'), ylabel('Energy')

figure, grid on, hold on
plot(t, Lsq)
xlabel('Time, seconds'), ylabel('Momentum')

figure, grid on, hold on
surface(xe, ye, ze)
title('Energy Ellipsoid')
axis('equal')
view([1,1,1])
xlabel('X'), ylabel('Y'), zlabel('Z')
colorbar

figure, grid on, hold on
surface(xm, ym, zm)
title('Momentum Ellipsoid')
axis('equal')
view([1,1,1])
xlabel('X'), ylabel('Y'), zlabel('Z')
colorbar

figure, grid on, hold on
surface(xp, yp, zp)
title('Polhode Ellipsoid')
axis('equal')
view([1,1,1])
colorbar

figure, grid on, hold on
%surface(xe-xm, ye-ym, ze-zm)
%surface(xm(1:11, :), ym(1:11, :,:), zm(1:11,:,:))
%surface(xm(11:21, :), ym(11:21, :,:), zm(11:21,:,:))
surface(xe, ye, ze)
surface(xm, ym, zm)
factor = 1;%15.6;
%plot3(factor*s(5,:), factor*s(6,:), factor*s(7,:), 'r','LineWidth', 5)
title('Energy and Momentum Ellipsoid')
axis('equal')
xlabel('X'), ylabel('Y'), zlabel('Z')
view([1,1,1])
%colorbar

% figure
% surface(omz)

figure, grid on, hold on
plot(a_wx, a_wy)
plot(wx, wy)
xlabel('\omega_x')
ylabel('\omega_y')
legend('Analytical', 'Numerical')
axis('equal')

figure, grid on, hold on
plot(t, error)
xlabel('Time, seconds')
ylabel('\Delta\omega, rad/s')
title('Difference in Angular Velocity between Numerical and Analytical')

figure, grid on, hold on
plot(t, a_wx)
%plot(t, a_wy)
plot(t, wx)
%plot(t, wy)
xlabel('Time, seconds')
ylabel('Angular Velocity, rad/s')