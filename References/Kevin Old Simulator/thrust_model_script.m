clear all
close all

t = linspace(0,8,1000);
T_des = 6*sin(t).^2;
dt = 0.1;
plot(t,T_des,'LineWidth',1)
T_real = zeros(size(T_des));

for i = 1:length(T_des)
    T_real(i+1) = thruster_model(T_des(i), dt, T_real(i));
end
hold on
T_real;
plot(t,T_real(2:end),'LineWidth',1)
title('Thruster Model')
legend('Commanded','Delivered')
xlabel('time (sec)')
ylabel('Thrust (N)')

%other test
T_des2 = (4.5/2)*square(3*t) + (4.5/2);
T_real2 = zeros(size(T_des2));
for i = 1:length(T_des)
    T_real2(i+1) = thruster_model(T_des2(i), dt, T_real2(i));
end
figure()
plot(t,T_des2,'LineWidth',1)
hold on
plot(t,T_real2(2:end),'LineWidth',1)
axis([0 max(t) -2 6])
title('Thruster Model')
legend('Commanded','Delivered')
xlabel('time (sec)')
ylabel('Thrust (N)')
