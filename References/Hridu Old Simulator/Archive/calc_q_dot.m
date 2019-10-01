function q_dot = calc_q_dot(q, w) 
%v1, not working so well
%     eta=q(1);
%     e1 =q(2);
%     e2 =q(3);
%     e3 =q(4);
%     e=[e1;e2;e3];   %euler axis
%     eta_dot=-1/2*e'*w;
%     e_dot= 1/2*(eta*eye(3)+Skew(e))*w; %euler angle
%     q_dot=[eta_dot; e_dot];

%v2, try DAmico's notes
wx = w(1); wy = w(2); wz = w(3);
OMEGA = [   0,    wz, -1*wy, wx;
        -1*wz,     0,    wx, wy;
           wy, -1*wx,     0, wz;
        -1*wx, -1*wy, -1*wz,  0];
q_dot = 0.5*OMEGA*q;
end