clc
clearvars
close all

GravityGradient_script
load SRP_Vectors.mat
load Feed_EKF.mat
M_totalgrav = M_sun + M_earth;
Model = [Qtn_model;Omega_model];

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, y_out] = ode45(@EulerEqsMomentsSRP_control, t, state0, options, M_totalgrav, tstep, c_sun,Rbus,Rsp,Nbus,Nsp,SAbus,SAsp,c_earth,Model);

A = zeros(3,3,n);
q = zeros(4,n);
A_bod2prin = [1     0       0;
                  0     0.9912  0.1325;
                  0     -0.1325  0.9912];
z_inertial = zeros(3,n);
for i = 1:n
    A(:,:,i) = q2A(y_out(i,7:10));
    z_inertial(:,i) = A(:,:,i)'*A_bod2prin*[0;0;1];
    z_desired = -c_earth(i,:)';
    q_des = zeros(4,1);
    q_des(1:3) = (z_inertial(:,i) + z_desired)/norm(z_inertial(:,i) + z_desired);
    A_error = A(:,:,i)*(q2A(q_des)*A(:,:,i))';
    q(:,i) = A2q(A_error);
end
[phi,theta,psi] = quat2EA312(q');

figure(1)
subplot(3,1,1)
plot(t_out,wrapToPi(phi))
ylabel('\phi (rad)')
subplot(3,1,2)
plot(t_out,wrapToPi(theta))
ylabel('\theta (rad)')
subplot(3,1,3)
plot(t_out,wrapToPi(psi))
xlabel('Time (s)')
ylabel('\psi (rad)')
% legend('\phi','\theta','\psi','location','best')
sgtitle('312 Euler Angles vs Time')

figure(2)
subplot(3,1,1)
plot(t_out,z_inertial(1,:))
ylabel('x')
subplot(3,1,2)
plot(t_out,z_inertial(2,:))
ylabel('y')
subplot(3,1,3)
plot(t_out,z_inertial(3,:))
ylabel('z')
xlabel('Time (s)')
% legend('\phi','\theta','\psi','location','best')
sgtitle('Z Axis in Inertial Frame vs Time')

save('dynamics_states.mat','t_out','y_out')