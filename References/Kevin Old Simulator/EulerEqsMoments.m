function [statedot] = EulerEqsMoments(t,state,M,tstep)
% state(1:3) are inertias
% state(4:6) are angular velocities; order is x, y, z
% state(7:10) are quaternion components
t_ind = int8(t/tstep) + 1;
% q1 = state(7);
% q2 = state(8);
% q3 = state(9);
% q4 = state(10);
% Q = [0 -q3 q2;
%         q3 0 -q1;
%        -q2 q1 0];
% q = [q1 q2 q3];
% A = (q4^2 - q'*q)*eye(3) + 2*(q*q') - 2*q4*Q;

statedot = zeros(10,1);
Ix = state(1);
Iy = state(2);
Iz = state(3);
wx = state(4);
wy = state(5);
wz = state(6);
statedot(4) = -(Iz - Iy)/Ix * wy * wz + M(t_ind,1)/Ix;
statedot(5) = -(Ix - Iz)/Iy * wz * wx + M(t_ind,2)/Iy;
statedot(6) = -(Iy - Ix)/Iz * wx * wy + M(t_ind,3)/Iz;
Omega = [0 wz -wy wx;
    -wz 0 wx wy;
    wy -wx 0 wz;
    -wx -wy -wz 0];
state(7:10) = state(7:10)/norm(state(7:10));
statedot(7:10) = 0.5*Omega*state(7:10);


end