function [mu, sig] = EKF_est(I, Model, c_sun, dt, t)
%function outputs estimatations of state and covariance matrix at each
%point in time. Inputs are the model values and the sun unit vector
%measurments. Angular velocities are assumed to come directly from the
%model. 

%moments of intertia
I = diag(I);
Ix = I(1);
Iy = I(2); 
Iz = I(3); 
K = [(Iz-Iy)/Ix, (Ix - Iz)/Iy, (Iy - Ix)/Iz];
S_sci = c_sun';
R = 0.001*eye(6);
Q_CV = 0.001*eye(7)*dt;

mu = [0 0 0 0 0 0]';
sig = 0.1*eye(7);

State = Model(:,1);
Qtn_model = Model(1:4);
Omega_model = Model(5:7);

Vt = mvnrnd(zeros(1,6),R,1)';
Wt = mvnrnd(zeros(1,7),Q_CV,1)';

f = [Qtn_model; Omega_model]; 

%calculate X
Xt = State;
Qtn = Xt(1:4);
Omega = Xt(5:7);
A = getAmatrix(t, Omega, K, Qtn);

%predict step, should happen before measurement step
[mu_pred, sig_pred] = predict(f, sig, A, Q_CV);

C = getHmatrix(Qtn, S_sci); 
Y = g(Xt,c_sun) + Vt; %determine measurement
Yhat = g(mu_pred,c_sun);

%update after measurement is taken 
[mu, sig] = update(mu_pred, sig_pred, C, R, Y, Yhat);
end

function Y = g(Xt,C_sun)
%Y = C*(Xt - mu_pred) + Vt; %determine measurement
q1 = Xt(1);
q2 = Xt(2);
q3 = Xt(3);
q4 = Xt(4);
wx = Xt(5);
wy = Xt(6);
wz = Xt(7);

Sx_sci = C_sun(1);
Sy_sci = C_sun(2);
Sz_sci = C_sun(3);

Sx = (q4^2 + q1^2 - q2^2 - q3^2)*Sx_sci + (2*q1*q2 + 2*q4*q3)*Sy_sci + (2*q1*q3 - 2*q4*q2)*Sz_sci;
Sy = (2*q1*q2 - 2*q4*q3)*Sx_sci + (q4^2 - q1^1 + q2^2 - q3^2)*Sy_sci + (2*q2*q3 + 2*q4*q1)*Sz_sci;
Sz = (2*q1*q3 + 2*q4*q2)*Sx_sci + (2*q2*q3 - 2*q4*q1)*Sy_sci + (q4^2 - q1^2 - q2^2 + q3^2)*Sz_sci;
Y = [Sx,Sy,Sz,wx,wy,wz]';
% Y = g + Vt; %determine measurement

end

function [mu_pred, sig_pred] = predict(f, sig, A, Q_CV)
    mu_pred = f;
    sig_pred = A*sig*A' + Q_CV;
end

function [mu, sig] = update(mu_pred, sig_pred, C, R, Y, g)    
    mu = mu_pred + sig_pred * C'*(C*sig_pred*C' + R)^-1 *(Y-g);
    sig = sig_pred - sig_pred*C'*(C*sig_pred*C' + R)^-1 *C*sig_pred;
end

function A = getAmatrix(t, W, K, Q)
%%Calculates A matrix from time, angular velocity, moments of inertia, and
%%quaternion values

%unpack quaternion
q1 = Q(1);
q2 = Q(2);
q3 = Q(3);
q4 = Q(4);

%unpack angular velocities
w_x = W(1);
w_y = W(2);
w_z = W(3);
w = norm([w_x w_y w_z]);

%unpack Ks (from Moments of Inertia)
k_x = K(1);
k_y = K(2);
k_z = K(3);

A = [cos((t*w)/2),  (w_z*sin((t*w)/2))/w, -(w_y*sin((t*w)/2))/w, (w_x*sin((t*w)/2))/w,...
    (q4*sin((t*w)/2))/w, -(q3*sin((t*w)/2))/w,  (q2*sin((t*w)/2))/w;
    -(w_z*sin((t*w)/2))/w, cos((t*w)/2),  (w_x*sin((t*w)/2))/w, (w_y*sin((t*w)/2))/w,...
    (q3*sin((t*w)/2))/w,  (q4*sin((t*w)/2))/w, -(q1*sin((t*w)/2))/w;
    (w_y*sin((t*w)/2))/w, -(w_x*sin((t*w)/2))/w, cos((t*w)/2), (w_z*sin((t*w)/2))/w,...
    -(q2*sin((t*w)/2))/w,  (q1*sin((t*w)/2))/w,  (q4*sin((t*w)/2))/w;
    -(w_x*sin((t*w)/2))/w, -(w_y*sin((t*w)/2))/w, -(w_z*sin((t*w)/2))/w, cos((t*w)/2),...
    -(q1*sin((t*w)/2))/w, -(q2*sin((t*w)/2))/w, -(q3*sin((t*w)/2))/w;
     0, 0, 0, 0, 1, k_x*t*w_x, k_x*t*w_y;
     0, 0, 0, 0,  k_y*t*w_z, 1, k_y*t*w_x;
     0, 0, 0, 0,   k_z*t*w_y, k_z*t*w_x, 1];
end

function H = getHmatrix(Qtn, S_sci)
%unpack quaternion
q1 = Qtn(1);
q2 = Qtn(2);
q3 = Qtn(3);
q4 = Qtn(4);

%unpack measurement
Sx_sci = S_sci(1);
Sy_sci = S_sci(2);
Sz_sci = S_sci(3);

H = [2*Sx_sci*q1 + 2*Sy_sci*q2 + 2*Sz_sci*q3, 2*Sy_sci*q1 - 2*Sx_sci*q2 - 2*Sz_sci*q4,...
    2*Sy_sci*q4 - 2*Sx_sci*q3 + 2*Sz_sci*q1, 2*Sx_sci*q4 + 2*Sy_sci*q3 - 2*Sz_sci*q2, 0, 0, 0;
    2*Sx_sci*q2 - Sy_sci + 2*Sz_sci*q4, 2*Sx_sci*q1 + 2*Sy_sci*q2 + 2*Sz_sci*q3,...
    2*Sz_sci*q2 - 2*Sy_sci*q3 - 2*Sx_sci*q4, 2*Sy_sci*q4 - 2*Sx_sci*q3 + 2*Sz_sci*q1, 0, 0, 0;
    2*Sx_sci*q3 - 2*Sy_sci*q4 - 2*Sz_sci*q1, 2*Sx_sci*q4 + 2*Sy_sci*q3 - 2*Sz_sci*q2,...
    2*Sx_sci*q1 + 2*Sy_sci*q2 + 2*Sz_sci*q3, 2*Sx_sci*q2 - 2*Sy_sci*q1 + 2*Sz_sci*q4, 0, 0, 0;
    0, 0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 0, 1];
end
