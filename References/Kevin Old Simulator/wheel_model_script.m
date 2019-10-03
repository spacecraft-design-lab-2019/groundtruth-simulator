clear all
close all

t = linspace(0,5,1000);
Tor_des = 5*sin(2*t);
dt = 0.1;
plot(t,Tor_des,'LineWidth',1)
Tor_real = zeros(size(Tor_des));
shift  = 0;
for i = 1:length(Tor_des)
    Tor_real(i+shift+1) = wheel_model(Tor_des(i), dt, Tor_real(i));
end
hold on
Tor_real;
plot(t,Tor_real(2:end-shift),'LineWidth',1)
title('Reaction/Momentum Wheel Model')
legend('Commanded','Delivered')
xlabel('time (sec)')
ylabel('Torque (N m)')
