function plotroe(tout, roe)
figure
subplot(3, 2, 1)
plot(tout, roe(:, 1), 'r')
title('a\delta a vs Time')
xlabel('Time [orbits]')
ylabel('a\delta a [km]') 
grid on

subplot(3, 2, 2)
plot(tout, roe(:, 2))
title('a\delta \lambda vs Time')
xlabel('Time [orbits]')
ylabel('a\delta \lambda [km]') 
grid on

subplot(3, 2, 3)
plot(tout, roe(:, 3))
title('a\delta e_x vs Time')
xlabel('Time [orbits]')
ylabel('a\delta e_x [km]') 
grid on

subplot(3, 2, 4)
plot(tout, roe(:, 4))
title('a\delta e_y vs Time')
xlabel('Time [orbits]')
ylabel('a\delta e_y [km]') 
grid on

subplot(3, 2, 5)
plot(tout, roe(:, 5))
title('a\delta i_x vs Time')
xlabel('Time [orbits]')
ylabel('a\delta i_x [km]') 
grid on

subplot(3, 2, 6)
plot(tout, roe(:, 6))
title('a\delta i_y vs Time')
xlabel('Time [orbits]')
ylabel('a\delta i_y [km]') 
grid on
end
