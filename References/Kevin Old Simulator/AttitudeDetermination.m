clc
clearvars -except psi_des
close all

FullDynamics_script

A_deter = zeros(3,3,n);
q_deter = zeros(4,n);
A_deter2 = zeros(3,3,n);
q_deter2 = zeros(4,n);
A_stat = zeros(3,3,n);

for i = 1:n
    % Deterministic Method
    v1 = v_sun(:,i);
    v2 = v_gyro(:,i);
    m1 = m_sun(:,i);
    pm = m1;
    m2 = m_gyro(:,i);
    qm = cross(m1,m2)/norm(cross(m1,m2));
    rm = cross(pm,qm);
    pv = v1;
    qv = cross(v1,v2)/norm(cross(v1,v2));
    rv = cross(pv,qv);
    
    M = [pm,qm,rm];
    V = [pv,qv,rv];
    
    A_deter(:,:,i) = M*inv(V);
    q_deter(:,i) = A2q(A_deter(:,:,i));
    
    % Ficticious measurement method
    m1_fict = (m1 + m2)/2;
    m1_fict = m1_fict/norm(m1_fict);
    v1_fict = (v1 + v2)/2;
    v1_fict = v1_fict/norm(v1_fict);
    v2_fict = (v1 - v2)/2;
    m2_fict = (m1 - m2)/2;
    pm_fict = m1_fict;
    qm_fict = cross(m1_fict,m2_fict)/norm(cross(m1_fict,m2_fict));
    rm_fict = cross(pm_fict,qm_fict);
    pv_fict = v1_fict;
    qv_fict = cross(v1_fict,v2_fict)/norm(cross(v1_fict,v2_fict));
    rv_fict = cross(pv_fict,qv_fict);
    
    M_fict = [pm_fict,qm_fict,rm_fict];
    V_fict = [pv_fict,qv_fict,rv_fict];
    
    A_deter2(:,:,i) = M_fict*inv(V_fict);
    q_deter2(:,i) = A2q(A_deter2(:,:,i));
end
[phi_deter,theta_deter,psi_deter] = quat2EA312(q_deter');
[phi_deter2,theta_deter2,psi_deter2] = quat2EA312(q_deter2');


q_stat = zeros(4,n);
for i = 1:n
    % Statistical Method
    v1 = v_sun(:,i);
    v2 = v_gyro(:,i);
    m1 = m_sun(:,i);
    m2 = m_gyro(:,i);
    w1 = 1;
    w2 = 1;
    W = [sqrt(w1)*m1, sqrt(w2)*m2];
    % U = [sqrt(w1)*v1, sqrt(w2)*v2];
    V = [v1, v2];
    B = W*V';
    S = B + B';
    Z = [B(2,3) - B(3,2), B(3,1) - B(1,3), B(1,2) - B(2,1)]';
    sigma = trace(B);
    K = [S - sigma*eye(3), Z; Z', sigma];
    [eigvecs,eigvals] = eig(K);
    [~,ind] = max(max(eigvals));
    q_stat(:,i) = eigvecs(:,ind);
    A_stat(:,:,i) = q2A(q(:,i));
end

[phi_stat,theta_stat,psi_stat] = quat2EA312(q_stat');

%% Desired Attitude Calculation

A_bod2prin = [1     0       0;
                  0     0.9912  0.1325;
                  0     -0.1325  0.9912];
q_des = zeros(4,n);
A_des = zeros(3,3,n);
z_test = zeros(3,n);
for i = 1:n
    z_inertial = A_stat(:,:,i)*A_bod2prin*[0;0;1];
    z_desired = -c_earth(i,:)';
    q_des(1:3,i) = (z_inertial + z_desired)/norm(z_inertial + z_desired);
    A_des(:,:,i) = q2A(q_des(:,i));
    z_test(:,i) = A_des(:,:,i)*z_inertial;
end

% save('A_des.mat','A_des')
[phi_des,theta_des,psi_des] = quat2EA312(q_des');


%% Plotting
figure(1)
subplot(3,1,1)
plot(t_out,wrapTo2Pi(phi))
hold on
plot(t_out,wrapTo2Pi(phi_deter))
ylabel('\phi (rad)')
legend('True','Est')
subplot(3,1,2)
plot(t_out,theta)
hold on
plot(t_out,theta_deter)
ylabel('\theta (rad)')
legend('True','Est')
subplot(3,1,3)
plot(t_out,psi)
hold on
plot(t_out,psi_deter)
xlabel('Time (s)')
ylabel('\psi (rad)')
legend('True','Est','location','best')
sgtitle('Deterministic Attitude Method Results with Ideal Sensors')

figure(2)
subplot(3,1,1)
plot(t_out,phi)
hold on
plot(t_out,wrapTo2Pi(phi_deter2))
ylabel('\phi (rad)')
legend('True','Est')
subplot(3,1,2)
plot(t_out,theta)
hold on
plot(t_out,theta_deter2)
ylabel('\theta (rad)')
legend('True','Est')
subplot(3,1,3)
plot(t_out,psi)
hold on
plot(t_out,psi_deter2 - 2*pi)
xlabel('Time (s)')
ylabel('\psi (rad)')
legend('True','Est','location','best')
sgtitle('Fictitious Deterministic Attitude Method Results')

figure(3)
subplot(3,1,1)
plot(t_out,phi)
hold on
plot(t_out,phi_stat)
ylabel('\phi (rad)')
legend('True','Est')
subplot(3,1,2)
plot(t_out,theta)
hold on
plot(t_out,theta_stat)
ylabel('\theta (rad)')
legend('True','Est')
subplot(3,1,3)
plot(t_out,psi)
hold on
plot(t_out,psi_stat-2*pi)
xlabel('Time (s)')
ylabel('\psi (rad)')
legend('True','Est','location','best')
sgtitle('Stastical Attitude Method Results with Ideal Sensors')

%% Error plots
figure(4)
subplot(3,1,1)
plot(t_out,wrapToPi(phi_deter-phi))
ylabel('\phi (rad)')
subplot(3,1,2)
plot(t_out,wrapToPi(theta_deter-theta))
ylabel('\theta (rad)')
subplot(3,1,3)
plot(t_out,wrapToPi(psi_deter-psi))
xlabel('Time (s)')
ylabel('\psi (rad)')
sgtitle('Deterministic Attitude Method Errors')

figure(5)
subplot(3,1,1)
plot(t_out,wrapToPi(phi_deter2-phi))
ylabel('\phi (rad)')
subplot(3,1,2)
plot(t_out,wrapToPi(theta_deter2-theta))
ylabel('\theta (rad)')
subplot(3,1,3)
plot(t_out,wrapToPi(psi_deter2-psi))
xlabel('Time (s)')
ylabel('\psi (rad)')
sgtitle('Fictious Deterministic Attitude Method Errors')

figure(6)
subplot(3,1,1)
plot(t_out,wrapToPi(phi_stat-phi))
ylabel('\phi (rad)')
subplot(3,1,2)
plot(t_out,wrapToPi(theta_stat-theta))
ylabel('\theta (rad)')
subplot(3,1,3)
plot(t_out,wrapToPi(psi_stat-psi))
xlabel('Time (s)')
ylabel('\psi (rad)')
sgtitle('Statistical Attitude Method Errors')

figure(7)
subplot(3,1,1)
plot(t_out,phi_des)
ylabel('\phi (rad)')
subplot(3,1,2)
plot(t_out,theta_des)
ylabel('\theta (rad)')
subplot(3,1,3)
plot(t_out,psi_des)
xlabel('Time (s)')
ylabel('\psi (rad)')
sgtitle('Attitude Control Error')
