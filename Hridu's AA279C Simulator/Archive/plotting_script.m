%plotting script for sim

t = (JD-JD(1))*(24*60*60);
%% position
figure, grid on, hold on
plot3(r_eci(1,:), r_eci(2,:), r_eci(3,:))
xlabel('X'),ylabel('Y'),zlabel('Z')
title('S/C Orbit, ECI')

figure, grid on, hold on
plot3(moon_eci(1,:), moon_eci(2,:), moon_eci(3,:))
xlabel('X'),ylabel('Y'),zlabel('Z')
title('Moon Orbit, ECI')

r_2_moon_eci = r_eci - moon_eci; 
figure, grid on, hold on
plot3(r_2_moon_eci(1,:), r_2_moon_eci(2,:), r_2_moon_eci(3,:))
xlabel('X'),ylabel('Y'),zlabel('Z')
title('S/C Orbit relative to Moon, ECI')

%% forces
figure
title('Accelerations, ECI')
subplot(4,1,1)
plot(t, a_eci_moon(1:3,:))
legend('x','y','z')
ylabel('Moon')
subplot(4,1,2)
plot(t, a_eci_earth(1:3,:))
ylabel('Earth')
subplot(4,1,3)
plot(t, a_eci_sun(1:3,:))
ylabel('Sun')
subplot(4,1,4)
plot(t, a_eci_srp(1:3,:))
ylabel('SRP')
xlabel('Time [s]')

