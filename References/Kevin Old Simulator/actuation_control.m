clear all
close all

load dynamics_states.mat

A = getA('xyz');
Ma = zeros(3,length(t_out));
wr = 6000 * 2*pi/60;
Ix = 184.9343576; % kg * m^2
Iy = 321.8711548; % kg * m^2
Iz = 423.2775247; % kg * m^2
I = diag([Ix Iy Iz]);

Ir = 12/wr;
Iw = diag([Ir Ir Ir]);
Ww = [wr wr wr]';
Lw = diag(Iw).*Ww;
Lwdot = [0 0 0]'; 

W = y_out(:,4:6);
W = W';
Wdot = diff(W,1,2);
Wdot = [[0,0,0]',Wdot];
Wdot(1,:) = lowpass(Wdot(1,:),0.001);
Wdot(2,:) = lowpass(Wdot(2,:),0.001);
Wdot(3,:) = lowpass(Wdot(3,:),0.001);

for i = 1:length(t_out)
    Ma(:,i) = torque_exchange(I, W(:,i), Wdot(:,i), A, Lw, Lwdot);
end

%% plotting
plot(t_out, Ma(1,:))
hold on
plot(t_out, Ma(2,:))
plot(t_out, Ma(3,:))
title('Actuator Modeled Torques')
xlabel('time (sec)')
ylabel('Torque (N-m)')