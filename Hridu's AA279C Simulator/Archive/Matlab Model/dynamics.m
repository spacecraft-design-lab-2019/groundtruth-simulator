%--------------------------------------------------------------------------
% SCRIPT: dynamics.m
% AUTHOR: Adam Zufall & Hridu Jain
% DATE:   April 23, 2019

% NOTES:
% Dynamics v2
%--------------------------------------------------------------------------

clc; clear all; close all;

%% Initialize
e_o = [0; 0; 0];            % initial Euler angles, rad
q_o = [0; 0; 0; 1];         % initial quaternion
%w_o = (pi/180)*[5;-3;7];    % initial rotation rate, rad/s
w_o = (pi/180)*[0.02; 5; -0.01];    % initial rotation rate, rad/s
s_o = [q_o; w_o];           % initial state, q w
s_o2= [e_o; w_o];

M = [0; 0; 0];

%% Moment of Inertia
[I, body_R_principal, ~] = moi_fxn();

Ix = I(1,1); Iy = I(2,2); Iz = I(3,3);
I_r = 3.44*10^-5;        %kg*m2,   moi of rotor
w_r = 500 * (2*pi/60); %rad/s, omega of rotor 5000 rpm
r_r = [0;1;0];          %   alignment of rotor
w_r_dot = 0;            %       constant rotor speed

%% Integrate Euler + Kinematic Equations
num_pts = 5*60*10;
t_span = linspace(0, 50*60, num_pts); % sim time = 5 min

% Quaternions
options = odeset('RelTol',1e-7);
[t, s] = ode45(@(t, s) rot_dyn(t, s, I, M, 'quaternions', w_r, w_r_dot, I_r, r_r), t_span, s_o, options);

% Euler Angles
[t_2, s_2] = ode45(@(t, s) rot_dyn(t, s, I, M, 'euler_angles', w_r, w_r_dot, I_r, r_r), t_span, s_o2, options);

%% Angular Momentum in Inertial Coordinates
run('orbit_driver.m')

q = [s(:,4), s(:,1:3)];
A = quat2dcm(q);

% Rotate w
w = s(:,end-2:end);
for i = 1:length(w)
    L_sc(i,:) = (A(:,:,i)' *(I * w(i,:)'))';
    L_rw(i,:) = (A(:,:,i)' *(I_r*w_r*r_r))';
    L(i, :) = ((A(:,:,i)' *(I * w(i,:)')) + A(:,:,i)' * (I_r*w_r*r_r))';
    w_inertial(i, :) = (A(:,:,i)' * w(i,:)')';
end


% %% Unit Vectors in 3D in Inertial Space
% 
% unit_inertial_x = [1;0;0];
% unit_inertial_y = [0;1;0];
% unit_inertial_z = [0;0;1];
% 
% for i = 1:length(t_out)
%     unit_orbital_x(i,:) = (rtn_R_eci(state_out(i,:), GM_Moon)'*unit_inertial_x)';
%     unit_orbital_y(i,:) = (rtn_R_eci(state_out(i,:), GM_Moon)'*unit_inertial_y)';
%     unit_orbital_z(i,:) = (rtn_R_eci(state_out(i,:), GM_Moon)'*unit_inertial_z)';
% end
% 
% for i = 1:length(A)
%     unit_principal_x(i,:) = (A(:,:,i)' * unit_inertial_x)';
%     unit_principal_y(i,:) = (A(:,:,i)' * unit_inertial_y)';
%     unit_principal_z(i,:) = (A(:,:,i)' * unit_inertial_z)';
%     
%     unit_body_x(i,:) = (body_R_principal * unit_principal_x(i,:)')';
%     unit_body_y(i,:) = (body_R_principal * unit_principal_y(i,:)')';
%     unit_body_z(i,:) = (body_R_principal * unit_principal_z(i,:)')';
% end
% 
% 
% %% Equilibrium Tests
% 
% % Given our particular orbit, we begin with RTN aligned with Moon Centered
% % Inertial Frame.
% 
% options2 = odeset('RelTol',1e-7, 'MaxStep',1e-2);
% 
% % Case 1: Assume Principal Axis aligned with Inertial Frame
% w_o3 = (pi/180)*[0;0;5];
% s_o3 = [e_o; w_o3];
% [t_3, s_3] = ode45(@(t, s) rot_dyn(t, s, I, M, 'euler_angles', w_r, w_r_dot, I_r, r_r), t_span, s_o3, options2);
% 
% 
% % Case 2: Set Initial Attitude to Match RTN Frame
% dcm_initial = rtn_R_eci(state_out(1, :), GM_Moon)'*eye(3); % Aligned with RTN, represented in inertial
% e_o4 = rotm2eul(dcm_initial)';
% s_o4 = [e_o4; w_o3];
% [t_4, s_4] = ode45(@(t, s) rot_dyn(t, s, I, M, 'euler_angles', w_r, w_r_dot, I_r, r_r), t_span, s_o4, options2);
% 
% testZer = s_4(:,1:3);
% for i = 1:length(s_4)
%     so = interp1(t_out, state_out, t_4(i));
%     A = eul2rotm(s_4(i, 1:3), 'XYZ');  %body2eci
%     R = rtn_R_eci(so, GM_Moon);  %eci2rtn
%     s_4(i, 1:3) = rotm2eul(A*R, 'XYZ')'; 
%     s_4(i, 4:6) = (R * s_4(i, 4:6)')';
%     
%     testOne(:,i) = rotm2eul(A, 'XYZ');
%     testTwo(:,i) = rotm2eul(R', 'XYZ');
%     testThr(:,i) = rotm2eul(A*R, 'XYZ');
% end

%% Plot

% Attitude - Quaternions
figure, hold on, grid on
plot(t, s(:,1:4))
legend('q1','q2','q3','q4')
title('Attitude - Quaternions')
xlabel('Time [s]')
ylabel('Quaternions')

% Attitude - Euler Angles
figure, hold on, grid on
plot(t_2, wrapTo180(rad2deg(s_2(:,1))))
plot(t_2, wrapTo180(rad2deg(s_2(:,2))))
plot(t_2, wrapTo360(rad2deg(s_2(:,3))))
legend('roll', 'pitch', 'yaw')
title('Attitude - Euler Angles')
xlabel('Time [s]')
ylabel('Euler Angles [deg]')

%Attitude - Q and Eul converted into Q
quatE = eul2quat(s_2(:,3:-1:1), 'ZYX');
figure, hold on, grid on
plot(t, s(:,1:4)-[quatE(:,2:4),quatE(:,1)])
xlabel('Time [s]'), ylabel('Error')
title('Difference between Quaternion Propagation and Euler Conversion')
legend('qx','qy','qz','q4')

%Attitude - Eul and Q convert to Eul
eulQ = quat2eul([s(:,4),s(:,1:3)], 'ZYX');
diff = wrapTo2Pi(s_2(:,1:3))-wrapTo2Pi(eulQ(:, 3:-1:1));
figure, hold on, grid on
plot(t_2, diff)
xlabel('Time [s]'), ylabel('Error [rad]')
title('Difference between Euler Angle Propagation and Quaternion Conversion')
legend('Roll','Pitch','Yaw')

% Angular Velocity Components
figure, hold on, grid on
plot(t, rad2deg(s(:,end-2:end)))
legend('w_x','w_y','w_z')
title('Angular Velocity Components')
xlabel('Time [s]')
ylabel('Angular Velocities [deg/s]')

% Angular Momentum Inertial Components
figure, hold on, grid on
subplot(3,1,1)
plot(t, L_sc(:,1:3))
title('Spacecraft Angular Momentum, Inertial Frame')
subplot(3,1,2)
plot(t, L_rw(:,1:3))
title('Rotor Angular Momentum, Inertial Frame')
ylabel('Angular Momentum Components [kg*m^2/sec]')
subplot(3,1,3)
plot(t, L(:,1:3))
legend('L_1', 'L_2', 'L_3')
title('Total Angular Momentum, Inertial Frame')
xlabel('Time [s]')

% 
% % Herpolode Angular Velocity Vector in Inertial
% figure, hold on, grid on,
% plot3(w_inertial(:,1),w_inertial(:,2),w_inertial(:,3))
% plot3([0, L(1,1)], [0, L(1,2)], [0, L(1,3)])
% axis([0 .1 0 .1 -.2 0])
% xlabel('x'), ylabel('y'), zlabel('z')
% x1 = w_inertial(1,1); y1 = w_inertial(1,2); z1 = w_inertial(1,3);
% w = null(L(1,:)); % Find two orthonormal vectors which are orthogonal to v
% [P,Q] = meshgrid(-50:50); % Provide a gridwork (you choose the size)
% X = x1+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
% Y = y1+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
% Z = z1+w(3,1)*P+w(3,2)*Q;
% surf(X,Y,Z)
% legend('Herpolode', 'Angular Momentum', 'Plane Normal to Angular Momentum')
% title('Herpolode Contained in Plane Perpendicular to Angular Momentum')
% 
% % Equilibrium Tests: Case 1 - Attitude and Velocities
% figure, hold on, grid on
% plot(t_3, wrapTo180(rad2deg(s_3(:,1:3))))
% legend('roll', 'pitch', 'yaw')
% title('Equilibrium Tests (Case 1): Attitude - Euler Angles')
% xlabel('Time [s]')
% ylabel('Euler Angles [deg]')
% 
% figure, hold on, grid on
% plot(t_3, rad2deg(s_3(:,end-2:end)))
% legend('w_x','w_y','w_z')
% title('Equilibrium Tests (Case 1): Angular Velocities')
% xlabel('Time [s]')
% ylabel('Angular Velocities [deg/s]')
% 
% % Equilibrium Tests: Case 2 - Attitude and Velocities
% figure, hold on, grid on
% plot(t_4, wrapTo180(rad2deg(s_4(:,1:3))))
% legend('roll', 'pitch', 'yaw')
% title('Equilibrium Tests (Case 2): Attitude - Euler Angles')
% xlabel('Time [s]')
% ylabel('Euler Angles [deg]')
% 
% figure, hold on, grid on
% plot(t_4, wrapTo180(wrapTo180(rad2deg(s_3(:,1:3)))-wrapTo180(rad2deg(s_4(:,1:3)))))
% xlabel('Time [s]'), ylabel('Difference [deg]')
% title('Difference in Euler Angles between Orbital and Inertial Frames')
% legend('Roll','Pitch','Yaw','Location','northwest')
% 
% figure, hold on, grid on
% plot(t_4, rad2deg(s_4(:,end-2:end)))
% legend('w_x','w_y','w_z')
% title('Equilibrium Tests (Case 2): Angular Velocities')
% xlabel('Time [s]')
% ylabel('Angular Velocities [deg/s]')
% 
% % Unit Vectors over Time
% figure, hold on, grid on
% pl(1) = plot3([0,1],[0,0],[0,0], '-b', 'LineWidth', 1.75);
% plot3([0,0],[0,1],[0,0], '-b', 'LineWidth', 1.75)
% plot3([0,0],[0,0],[0,1], '-b', 'LineWidth', 1.75)
% for tsteps = [1, round(.25*length(t_out)), round(.5*length(t_out)), round(.75*length(t_out)), length(t_out)]
%     pl(2) = plot3([0,unit_orbital_x(tsteps,1)],[0,unit_orbital_x(tsteps,2)],[0,unit_orbital_x(tsteps,3)], '--r');
%     pl(3) = plot3([0,unit_orbital_y(tsteps,1)],[0,unit_orbital_y(tsteps,2)],[0,unit_orbital_y(tsteps,3)], '--g');
%     pl(4) = plot3([0,unit_orbital_z(tsteps,1)],[0,unit_orbital_z(tsteps,2)],[0,unit_orbital_z(tsteps,3)], '--c');
% end
% axis([-1 1 -1 1 -1 1])
% axis('equal')
% xlabel('x'), ylabel('y'), zlabel('z')
% legend([pl(1:4)], 'Inertial', 'Orbital Frame X', 'Orbital Frame Y', 'Orbital Frame Z')
% title('Unit Vectors of Orbital Frame over Time')
% 
% figure, hold on, grid on
% pl(1) = plot3([0,1],[0,0],[0,0], '-b', 'LineWidth', 1.75);
% plot3([0,0],[0,1],[0,0], '-b', 'LineWidth', 1.75)
% plot3([0,0],[0,0],[0,1], '-b', 'LineWidth', 1.75)
% for tsteps = [1, 4, 8, 12]%round(.025*length(t)), round(.05*length(t)), round(.075*length(t)), length(t)]
%     pl(2) = plot3([0,unit_principal_x(tsteps,1)],[0,unit_principal_x(tsteps,2)],[0,unit_principal_x(tsteps,3)], '--r');
%     pl(3) = plot3([0,unit_principal_y(tsteps,1)],[0,unit_principal_y(tsteps,2)],[0,unit_principal_y(tsteps,3)], '--g');
%     pl(4) = plot3([0,unit_principal_z(tsteps,1)],[0,unit_principal_z(tsteps,2)],[0,unit_principal_z(tsteps,3)], '--c');
% end
% axis([-1 1 -1 1 -1 1])
% axis('equal')
% xlabel('x'), ylabel('y'), zlabel('z')
% legend([pl(1:4)], 'Inertial', 'Principal Frame X', 'Principal Frame Y', 'Principal Frame Z')
% title('Unit Vectors of Principal Frame over Time')
% 
% figure, hold on, grid on
% pl(1) = plot3([0,1],[0,0],[0,0], '-b', 'LineWidth', 1.75);
% plot3([0,0],[0,1],[0,0], '-b', 'LineWidth', 1.75)
% plot3([0,0],[0,0],[0,1], '-b', 'LineWidth', 1.75)
% for tsteps = [1, 4, 8, 12]%round(.025*length(t)), round(.05*length(t)), round(.075*length(t)), length(t)]
%     pl(2) = plot3([0,unit_body_x(tsteps,1)],[0,unit_body_x(tsteps,2)],[0,unit_body_x(tsteps,3)], '--r');
%     pl(3) = plot3([0,unit_body_y(tsteps,1)],[0,unit_body_y(tsteps,2)],[0,unit_body_y(tsteps,3)], '--g');
%     pl(4) = plot3([0,unit_body_z(tsteps,1)],[0,unit_body_z(tsteps,2)],[0,unit_body_z(tsteps,3)], '--c');
% end
% axis([-1 1 -1 1 -1 1])
% axis('equal')
% xlabel('x'), ylabel('y'), zlabel('z')
% legend([pl(1:4)], 'Inertial', 'Body Frame X', 'Body Frame Y', 'Body Frame Z')
% title('Unit Vectors of Body Frame over Time')