clc
clearvars
close all

GravityGradient_script
load SRP_Vectors.mat

M_totalgrav = M_sun + M_earth;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
[t_out, y_out] = ode45(@EulerEqsMomentsSRP, t, state0, options, M_totalgrav, tstep, c_sun,Rbus,Rsp,Nbus,Nsp,SAbus,SAsp);


%% Plotting

q1t = y_out(:,7);
q2t = y_out(:,8);
q3t = y_out(:,9);
q4t = y_out(:,10);

% phi = unwrap(atan2(2*q1t.*q3t - 2*q4t.*q2t,2*q2t.*q3t+ 2*q4t.*q1t));
% theta = acos(q3t.^2 - q2t.^2 - q1t.^2 + q4t.^2);
% psi = unwrap(atan2(2*q1t.*q3t+2*q4t.*q2t,-2*q2t.*q3t+ 2*q4t.*q1t));
% 
% A_prin2bod = [1     0       0;
%                   0     0.9912  -0.1325;
%                   0     0.1325  0.9912];
A = zeros(3,3,length(t));
m_sun_val = zeros(3,length(t));
m_gyro_val = zeros(3,length(t));
q = zeros(4,length(t));
for i = 1:length(t)
    q(:,i) = [q1t(i), q2t(i), q3t(i), q4t(i)];
    A(:,:,i) = q2A(q(:,i));
    m_sun_val(:,i) = A(:,:,i)*[cxsun(i); cysun(i); czsun(i)] + 0*(mvnrnd(zeros(3,1),1e-7*eye(3))' + 1e-4*[1; 1; 1]);
    m_sun_val(:,i) = m_sun_val(:,i)/norm(m_sun_val(:,i));
    m_gyro_val(:,i) = A(:,:,i)*(y_out(i,4:6)') + 0*( mvnrnd(zeros(3,1),1e-7*eye(3))' + 1e-4*[1; 1; 1]);
    m_gyro_val(:,i) = m_gyro_val(:,i)/norm(m_gyro_val(:,i));
end
m_gyro = m_gyro_val;
m_sun = m_sun_val;
v_sun = [cxsun, cysun, czsun]';
v_gyro = y_out(:,4:6)';

[phi,theta,psi] = quat2EA312(q');
% m_sun.time = [];
% m_sun.signals.values = m_sun_val';
% m_sun.signals.dimensions = 3;
% m_gyro.time = [];
% m_gyro.signals.values = m_gyro_val';
% m_gyro.signals.dimensions = 3;
% v_sun.time = [];
% v_sun.signals.values = [cxsun, cysun, czsun];
% v_sun.signals.dimensions = 3;
% v_gyro.time = [];
% v_gyro.signals.values = y_out(:,4:6);
% v_gyro.signals.dimensions = 3;
% 
% m1 = m_sun.signals.values;
% m2 = m_gyro.signals.values;
% ms = [m1, m2];
% v1 = v_sun.signals.values;
% v2 = v_gyro.signals.values;
% vs = [v1, v2];
% ts = t;
% As = A;
% save('myData2.mat', 'ms','vs','ts','As')
% 
% %%
% clc
% close all
% 
% figure(1)
% subplot(3,1,1)
% plot(t(1:200),reshape(A(1,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_deter(1,1,1:200),[1,200]))
% subplot(3,1,2)
% plot(t(1:200),reshape(A(2,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_deter(2,1,1:200),[1,200]))
% subplot(3,1,3)
% plot(t(1:200),reshape(A(3,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_deter(3,1,1:200),[1,200]))
% 
% figure(2)
% subplot(3,1,1)
% plot(t(1:200),reshape(A(1,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(1,1,1:200),[1,200]))
% subplot(3,1,2)
% plot(t(1:200),reshape(A(2,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(2,1,1:200),[1,200]))
% subplot(3,1,3)
% plot(t(1:200),reshape(A(3,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(3,1,1:200),[1,200]))
% 
% figure(3)
% subplot(3,1,1)
% plot(t(1:200),reshape(A(1,2,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(1,2,1:200),[1,200]))
% subplot(3,1,2)
% plot(t(1:200),reshape(A(2,2,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(2,2,1:200),[1,200]))
% subplot(3,1,3)
% plot(t(1:200),reshape(A(3,2,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(3,2,1:200),[1,200]))
% 
% figure(4)
% subplot(3,1,1)
% plot(t(1:200),reshape(A(1,3,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(1,3,1:200),[1,200]))
% subplot(3,1,2)
% plot(t(1:200),reshape(A(2,3,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(2,3,1:200),[1,200]))
% subplot(3,1,3)
% plot(t(1:200),reshape(A(3,3,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_stat(3,3,1:200),[1,200]))
% 
% figure(5)
% subplot(4,1,1)
% plot(t,q(1,:))
% hold on
% plot(t,q_out(:,1))
% subplot(4,1,2)
% plot(t,q(2,:))
% hold on
% plot(t,q_out(:,2))
% subplot(4,1,3)
% plot(t,q(3,:))
% hold on
% plot(t,q_out(:,3))
% subplot(4,1,4)
% plot(t,q(4,:))
% hold on
% plot(t,q_out(:,4))
% 
% [phi, theta, psi] = quat2EA(q');
% [phi2, theta2, psi2] = quat2EA(q_out);
% 
% figure(6)
% subplot(3,1,1)
% plot(t_out,phi)
% hold on
% plot(t_out,phi2)
% subplot(3,1,2)
% plot(t_out,theta)
% hold on
% plot(t_out,theta2)
% subplot(3,1,3)
% plot(t_out,psi)
% hold on
% plot(t_out,psi2)
% xlabel('Time (s)')
% ylabel('Angle (rad)')
% % legend('\phi','\theta','\psi','location','best')
% sgtitle('313 Euler Angles vs Time')
% 
% %% find best q values
% clc
% 
% % q_norm = zeros(1,4);
% % q_best = zeros(4,length(t));
% % for i = 1:length(t)
% %     for j = 1:4
% %         q_norm(j) = norm(eigvecs(:,j,i) - A2q(A_deter(:,:,i)));
% %     end
% %     [~,min_ind] = min(q_norm);
% %     q_best(:,i) = eigvecs(:,min_ind,i);
% % end
% A_norm = zeros(1,4);
% A_best = zeros(3,3,length(t));
% for i = 1:length(t)
%     for j = 1:4
%         Atest = q2A(eigvecs(:,j,i));
%         A_norm(j) = sum(sum(Atest-A_deter(:,:,i)));
%     end
%     [~,min_ind] = min(A_norm);
%     A_best(:,:,i) = q2A(eigvecs(:,min_ind,i));
% end
% % 
% % figure(7)
% % subplot(4,1,1)
% % plot(t,q(1,:))
% % hold on
% % plot(t,q_best(1,:))
% % subplot(4,1,2)
% % plot(t,q(2,:))
% % hold on
% % plot(t,q_best(2,:))
% % subplot(4,1,3)
% % plot(t,q(3,:))
% % hold on
% % plot(t,q_best(3,:))
% % subplot(4,1,4)
% % plot(t,q(4,:))
% % hold on
% % plot(t,q_best(4,:))
% 
% figure(8)
% subplot(3,1,1)
% plot(t(1:200),reshape(A(1,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_best(1,1,1:200),[1,200]))
% subplot(3,1,2)
% plot(t(1:200),reshape(A(2,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_best(2,1,1:200),[1,200]))
% subplot(3,1,3)
% plot(t(1:200),reshape(A(3,1,1:200),[1,200]))
% hold on
% plot(t(1:200),reshape(A_best(3,1,1:200),[1,200]))