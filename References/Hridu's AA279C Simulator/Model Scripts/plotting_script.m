%plotting script for sim

t = (JD-JD(1))*(24*60*60);
ang = (pi/180) * (SS.constants.Obliquity_Earth - SS.constants.Obliquity_Moon);
R_eci2mci = [1,0,0;
             0,cos(ang),-1*sin(ang);
             0,sin(ang),   cos(ang)];
%% position
figure, grid on, hold on
plot3(r_eci(:,1), r_eci(:,2), r_eci(:,3))
xlabel('X'),ylabel('Y'),zlabel('Z')
title('S/C Orbit, ECI')
axis('equal')

figure, grid on, hold on
plot3(moon_eci(:,1), moon_eci(:,2), moon_eci(:,3))
xlabel('X'),ylabel('Y'),zlabel('Z')
title('Moon Orbit, ECI')
axis('equal')

r_mci = (R_eci2mci*(r_eci - moon_eci)')'; 
[mx,my,mz] = ellipsoid(0,0,0,1738100,1738100,1736000,50);
figure, grid on, hold on
plot3(r_mci(:,1), r_mci(:,2), r_mci(:,3), 'LineWidth',1)
surf(mx,my,mz);
props.FaceColor= [0.6 0.6 0.6];
% s.FaceColor.red = 0.6;
% s.FaceColor.blue = 0.6;
% s.FaceColor.green = 0.6;
xlabel('X'),ylabel('Y'),zlabel('Z')
title('S/C Orbit relative to Moon, MCI')
axis('equal')

for i=1:length(t)
    temp(i) = norm(r_mci(i,:));
    mechEnergy_eci(i) = 0.5*norm(v_eci(i,:))^2 - SS.constants.GM_Earth/norm(r_eci(i,:));
    v_mci = R_eci2mci * (v_eci(i,:)' - v_moon_eci(i,:)');
    mechEnergy_mci(i) = 0.5*norm(v_mci)^2 - SS.constants.GM_Moon /norm(r_mci(i,:));
end
figure, grid on, hold on
plot(t/60, temp)
xlabel('Time [min]'), ylabel('|R_{MCI}|')
title('Distance to Moon')

figure, grid on, hold on
plot(t/60, mechEnergy_eci)
xlabel('Time [min]'), ylabel('Energy')
title('Specific Mechanical Energy, ECI')

figure, grid on, hold on
plot(t(2:end)/60, mechEnergy_mci(2:end))
xlabel('Time [min]'), ylabel('Energy')
title('Specific Mechanical Energy, MCI')

%% forces
figure
subplot(4,1,1)
plot(t/60, a_eci_moon(:,1:3))
legend('x','y','z')
title('Accelerations, ECI, m/s^2')
ylabel('Moon')
subplot(4,1,2)
plot(t/60, a_eci_earth(:,1:3))
ylabel('Earth')
subplot(4,1,3)
plot(t/60, a_eci_sun(:,1:3))
ylabel('Sun')
subplot(4,1,4)
plot(t/60, a_eci_srp(:,1:3))
ylabel('SRP')
xlabel('Time [min]')

%% torques
figure
subplot(4,1,1)
plot(t/60, tor_gg_earth(:,1:3))
legend('x','y','z')
title('Torques, Body, N*m')
ylabel('GG Earth')
subplot(4,1,2)
plot(t/60, tor_gg_moon(:,1:3))
ylabel('Moon')
subplot(4,1,3)
plot(t/60, tor_srp(:,1:3))
ylabel('SRP')
subplot(4,1,4)
plot(t/60, tor_tot_body(:,1:3))
ylabel('Total')
xlabel('Time [min]')

%% attitude

for i=1:length(t)  
    rol(i) = (180/pi)*wrapTo2Pi(e_out(i,1));
    pit(i) = (180/pi)*wrapTo2Pi(e_out(i,2));
    yaw(i) = (180/pi)*wrapTo2Pi(e_out(i,3));
    if(rol(i) > 180)
        rol(i) = rol(i)-360;
    end
    if(pit(i) > 180)
        pit(i) = pit(i)-360;
    end
    if(yaw(i) > 180)
        yaw(i) = yaw(i)-360;
    end
end
figure, grid on, hold on
%plot(t, (180*pi)*(e_out(:,1:3)))
plot(t/60, rol,'.','MarkerSize',3)
plot(t/60, pit,'.','MarkerSize',3)
plot(t/60, yaw,'.','MarkerSize',3)
xlabel('Time [min]')
ylabel('Angle [deg]')
title('Euler Angle Time History, ECI to Body')
legend('Roll','Pitch','Yaw')

% e_true = squeeze(e_true);
% e_true = e_true';
% rol = (180/pi)*(e_true(:,1));
% pit = (180/pi)*(e_true(:,2));
% for i=1:length(t)    
%     yaw(i) = (180/pi)*wrapTo2Pi(e_true(i,3));
%     if(yaw(i) > 180)
%         yaw(i) = yaw(i)-360;
%     end
% end
% figure, grid on, hold on
% %plot(t, (180*pi)*(e_out(:,1:3)))
% plot(t/60, rol,'.','MarkerSize',3)
% plot(t/60, pit,'.','MarkerSize',3)
% plot(t/60, yaw,'.','MarkerSize',3)
% xlabel('Time [min]')
% ylabel('Angle [deg]')
% title('Euler Angle Time History, MCI to Body')
% legend('Roll','Pitch','Yaw')

%% omega
figure
plot(t/60, w_MCI(:,1:3))
title('Omega MCI')

figure
plot(t/60, w_body(:,1:3))
title('Omega Body')

%% antennae
figure, grid on, hold on
plot(t/60, ant_axis_MCI(:,1:3))
title('Ant')

for i=1:length(t)
    dotp(i) = dot([0;0;1], ant_axis_MCI(i,1:3));
end
ang_error = acosd(dotp);
figure, grid on, hold on
plot(t/60, real(ang_error))
title('Antenae Axis Pointing Error')
ylabel('Angle [deg]')
xlabel('Time [min]')


%% HW 6 Plots
% figure, grid on, hold on
% plot(t/60, q_out(:,1:4) - q_est_det(:,1:4))
% xlabel('Time [min]'), ylabel('Error')
% title('Attitude Error, With Sensor Noise')
% legend('q1','q2','q3','q4')

% A = quat2dcm([q_out(:,4),q_out(:,1:3)]);
% B = quat2dcm([q_est_int(:,4),q_est_int(:,1:3)]);
% for(i=1:length(t))
%     C(:,:,i) = A(:,:,i)'*B(:,:,i);
% end
% eul_diff = (180/pi)*rotm2eul(C, 'XYZ');
% figure, grid on, hold on
% plot(t/60, eul_diff)
% xlabel('Time [min]'), ylabel('Error')
% title('Attitude Error, Integrating Rates, With Sensor Noise')
% legend('Roll','Pitch','Yaw')
% 
% A = quat2dcm([att_estimate(:,4),att_estimate(:,1:3)]);
% B = quat2dcm([att_desired(:,4),att_desired(:,1:3)]);
% for(i=1:length(t))
%     C(:,:,i) = A(:,:,i)'*B(:,:,i);
% end
% eul_diff = (180/pi)*rotm2eul(C, 'XYZ');
% figure, grid on, hold on
% plot(t/60, eul_diff)
% xlabel('Time [min]'), ylabel('Error')
% title('Difference Between Desired and Estimated Attitude, with Sensor Noise')
% legend('Roll','Pitch','Yaw')

%% HW 7 Plots
figure, grid on, hold on
plot(t/60, q_out(:,1:4) - q_EKF(:,1:4))
xlabel('Time [min]'), ylabel('Error')
title('True Attitude Estimation Error, from EKF')
legend('q1','q2','q3','q4')

A = quat2dcm([q_out(:,4),q_out(:,1:3)]);
B = quat2dcm([q_EKF(:,4),q_EKF(:,1:3)]);
for(i=1:length(t))
    C(:,:,i) = A(:,:,i)'*B(:,:,i);
end
% eul_diff = (180/pi)*rotm2eul(C, 'XYZ');
% figure, grid on, hold on
% plot(t/60, eul_diff)
% xlabel('Time [min]'), ylabel('Error [deg]')
% title('True Attitude Estimation Error, from EKF')
% legend('Roll','Pitch','Yaw')

figure, grid on, hold on
plot(t/60, omega_out(:,1:3) - w_EKF(:,1:3))
xlabel('Time [min]'), ylabel('Error [rad/s]')
title('True Angular Velocity Estimation Error, from EKF')
legend('\omega_x','\omega_y','\omega_z')

figure, grid on, hold on
plot(t/60, squeeze(sigma_EKF(1,1,:)))
plot(t/60, squeeze(sigma_EKF(2,2,:)))
plot(t/60, squeeze(sigma_EKF(3,3,:)))
plot(t/60, squeeze(sigma_EKF(4,4,:)))
xlabel('Time [min]'), ylabel('Error')
title('Estimated Attitude Estimation Error, from EKF')
legend('q1','q2','q3','q4')

% figure, grid on, hold on
% plot(t/60, y_prefit_EKF(:,1:3),'.')
% plot(t/60, y_postfit_EKF(:,1:3),'.')
% xlabel('Time [min]'), ylabel('Residuals (Magnitude)')
% legend('Pre_x','Pre_y','Pre_z','Post_x','Post_y','Post_z')
% title('\Delta Sun Measurements')

% figure, grid on, hold on
% plot(t/60, y_prefit_EKF(:,4:6),'.')
% plot(t/60, y_postfit_EKF(:,4:6),'.')
% xlabel('Time [min]'), ylabel('Residuals [rad/s]')
% legend('Pre_x','Pre_y','Pre_z','Post_x','Post_y','Post_z')
% title('\Delta Angular Rates')

% figure, grid on, hold on
% plot(t/60, vecnorm(y_prefit_EKF(:,1:3),2,2),'.')
% plot(t/60, vecnorm(y_postfit_EKF(:,1:3),2,2),'.')
% xlabel('Time [min]'), ylabel('Residuals (Magnitude)')
% title('Sun Direction Residuals')
% legend('Pre','Post')

% figure, grid on, hold on
% plot(t/60, vecnorm(y_prefit_EKF(:,4:6),2,2),'.')
% plot(t/60, vecnorm(y_postfit_EKF(:,4:6),2,2),'.')
% xlabel('Time [min]'), ylabel('Residuals (Magnitude)')
% title('Angular Rate Residuals')
% legend('Pre','Post')

%% HW 8 Plots
figure, grid on, hold on
plot(t/60, act_tor_body(:,1:3))
xlabel('Time [min]')
ylabel('Torque [Nm]')
title('Actuator Torques in Body Frame')
legend('X','Y','Z')

figure, grid on, hold on
plot(t/60, w_EKF(:,1:3))
xlabel('Time [min]')
ylabel('Angular Rate [rad/s]')
title('Estimated Angular Rate, Body Frame')
legend('X','Y','Z')

% figure, grid on, hold on
% plot(t/60, act_tor_body - tor_act_est)
% xlabel('Time [min]')
% ylabel('Error [Nm]')
% title('Actuator Deviation from Ideal Torque, Body Frame')
% legend('X','Y','Z')

figure, grid on, hold on
plot(t/60, omega_out(:,1:3))
xlabel('Time [min]')
ylabel('Rate [rad/s]')
title('Actual Angular Velocity, Body Frame')

figure, grid on, hold on
plot(t/60, omega_out(:,1:3)-gyro_meas(:,1:3))
xlabel('Time [min]')
ylabel('Rate [rad/s]')
title('Gyro Measurements Angular Velocity Error, Body Frame')

%% final plots 
% *Europe plays faintly in the background*

mode_adj_exp_tor = zeros(length(t),3);
mode_adj_ang_vel = zeros(length(t),3);
for i=1:length(t)
    if(mode_out(i) == 1)
        mode_adj_exp_tor(i,:) = tor_act_est(i,:);
        mode_adj_ang_vel(i,:) = w_EKF(i,:);
    else
        mode_adj_exp_tor(i,:) = tor_act_est1(i,:);
        mode_adj_ang_vel(i,:) = w_EKF1(i,:);
    end
end

figure, grid on, hold on
plot(t/60, act_tor_body - mode_adj_exp_tor)
xlabel('Time [min]')
ylabel('Error [Nm]')
title('Actuator Deviation from Ideal Torque, Body Frame')
legend('X','Y','Z')

figure, grid on, hold on
plot(t/60, omega_out - mode_adj_ang_vel)
xlabel('Time [min]')
ylabel('Error [rad/s]')
title('Actual - Estimated Angular Velocity, Body Frame')
legend('X','Y','Z')


%% ellipsoids -- we all hate them but we need them

I = SC.mass_properties.inertia;

figure, grid on, hold on
plot3(w_body(:,1), w_body(:,2), w_body(:,3))
xlabel('X'), ylabel('Y'), zlabel('Z')
title('Angular Velocity Profile')

wx = w_body(:,1);
wy = w_body(:,2);
wz = w_body(:,3);

Ix = I(1,1); Iy = I(2,2); Iz = I(3,3);

for i=1:length(t)
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

% % rxp = ((I(2,2)-Lsq(1)*twoT(1))*(I(3,3)-Lsq(1)*twoT(1)))^-1;
% % ryp = ((I(1,1)-Lsq(1)*twoT(1))*(I(3,3)-Lsq(1)*twoT(1)))^-1;
% % rzp = ((I(2,2)-Lsq(1)*twoT(1))*(I(1,1)-Lsq(1)*twoT(1)))^-1;
% % [xp, yp, zp] = ellipsoid(0,0,0, rxp, ryp, rzp, num);

% % omx = linspace(min(wx), max(wx), num);
% % omy = linspace(min(wy), max(wy), num);
% % 
% % L2 = mean(Lsq);
% % T2 = mean(twoT);
% % for a = 1:num
% %     for b = 1:num
% %         p1 = -1*(Ix-L2/T2)*Ix*omx(a)*omx(a);
% %         p2 = -1*(Iy-L2/T2)*Iy*omy(b)*omy(b);
% %         p3 = (Iz-L2/T2)*Iz;
% %         omz(a,b) = sqrt((p1+p2)/p3);
% %     end
% % end

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

figure, grid on, hold on
surface(xm, ym, zm)
title('Momentum Ellipsoid')
axis('equal')
view([1,1,1])
xlabel('X'), ylabel('Y'), zlabel('Z')


figure, grid on, hold on
surface(xe, ye, ze)
surface(xm, ym, zm)
plot3(wx, wy, wz, 'r', 'LineWidth', 9)
title('Energy and Momentum Ellipsoid [kg*m^2/s^2]')
axis('equal')
xlabel('X'), ylabel('Y'), zlabel('Z')
view([1,1,1])



