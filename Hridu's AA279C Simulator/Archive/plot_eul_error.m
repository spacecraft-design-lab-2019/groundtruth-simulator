%% attitude
t = 0:0.05:0.05*(length(rol)-1)
rol = (180/pi)*(e_out_w_torques(:,1)-e_out_wo_torques(:,1));
pit = (180/pi)*(e_out_w_torques(:,2)-e_out_wo_torques(:,2));
yaw = (180/pi)*(e_out_w_torques(:,3)-e_out_wo_torques(:,3));
% for i=1:length(t)    
%     yaw(i) = (180/pi)*wrapTo2Pi(e_out(i,3));
%     if(yaw(i) > 180)
%         yaw(i) = yaw(i)-360;
%     end
% end
figure, grid on, hold on
%plot(t, (180*pi)*(e_out(:,1:3)))
plot(t/60, rol,'.','MarkerSize',3)
plot(t/60, pit,'.','MarkerSize',3)
plot(t/60, yaw,'.','MarkerSize',3)
xlabel('Time [min]')
ylabel('Angle [deg]')
title('Euler Angle Difference, ECI to Body')
legend('Roll','Pitch','Yaw')