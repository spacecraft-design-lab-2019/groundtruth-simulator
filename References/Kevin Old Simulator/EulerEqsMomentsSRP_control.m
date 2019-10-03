function [statedot] = EulerEqsMomentsSRP_control(t,state,M,tstep,c_sun,Rbus,Rsp,Nbus,Nsp,SAbus,SAsp,c_earth,Model)
% state(1:3) are inertias
% state(4:6) are angular velocities; order is x, y, z
% state(7:10) are quaternion components
t_ind = int8(t/tstep) + 1;
q1 = state(7);
q2 = state(8);
q3 = state(9);
q4 = state(10);
Q = [0 -q3 q2;
        q3 0 -q1;
       -q2 q1 0];
A = [q1^2 - q2^2 - q3^2 + q4^2, 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4);
    2*(q1*q2 - q3*q4), -q1^2 + q2^2 - q3^2 + q4^2, 2*(q2*q3 + q1*q4);
    2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), -q1^2 - q2^2 + q3^2 + q4^2];
A_prin2bod = [1     0       0;
              0     0.9912  -0.1325;
              0     0.1325  0.9912];
s = -A_prin2bod*A*(c_sun(t_ind,:)');
s = s/norm(s);

CS_bus = 0.2;
CD_bus = 0.7;
CS_sp = 0.35;
P = 1358/299792; % momentum flux = F/c

M_SRP = zeros(3,1);
for i = 1:8
    lit = Nbus(:,i)'*s;
    if lit > 0
        F = -P*SAbus(i)*lit*((1-CS_bus)*s + 2*(CS_bus*lit + 1/3 * CD_bus)*Nbus(:,i));
        r = Rbus(:,i);
        M_SRP = M_SRP + cross(r,F);
    end
end
for j = 1:12
    lit = Nsp(:,j)'*s;
    if lit > 0
        F = -P*SAsp(j)*lit*((1-CS_sp)*s + 2*CS_sp*lit*Nsp(:,j));
        r = Rsp(:,j);
        M_SRP = M_SRP + cross(r,F);
    end
end
Ix = state(1);
Iy = state(2);
Iz = state(3);
wx = state(4);
wy = state(5);
wz = state(6);
I = diag([Ix Iy Iz]);

[Mu, ~] = EKF_est(I, Model(:,t_ind), c_sun(t_ind,:), tstep, t);
W = Mu(5:7);

M_SRP_prin = A_prin2bod' * M_SRP;
z_inertial = A'*A_prin2bod'*[0;0;1];
z_desired = -c_earth(t_ind,:)';
q_des = zeros(4,1);
q_des(1:3) = (z_inertial + z_desired)/norm(z_inertial + z_desired);
A_des = q2A(q_des)*A;
A_E = A*A_des';


M_c = control_torque(I, W, M(t_ind,:)', A_E)';
M_tot = M(t_ind,:)' + M_SRP_prin + M_c;
statedot = zeros(10,1);




statedot(4) = -(Iz - Iy)/Ix * wy * wz + M_tot(1)/Ix;
statedot(5) = -(Ix - Iz)/Iy * wz * wx + M_tot(2)/Iy;
statedot(6) = -(Iy - Ix)/Iz * wx * wy + M_tot(3)/Iz;
Omega = [0 wz -wy wx;
    -wz 0 wx wy;
    wy -wx 0 wz;
    -wx -wy -wz 0];
state(7:10) = state(7:10)/norm(state(7:10));
statedot(7:10) = 0.5*Omega*state(7:10);

%save('M_SRP.mat','M_SRP_prin','-append')

end