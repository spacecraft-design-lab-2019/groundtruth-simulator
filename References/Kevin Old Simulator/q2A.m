function A = q2A(q)
q1 = q(1);
q2 = q(2);
q3 = q(3);
q4 = q(4);
% Q = [0 -q3 q2;
%         q3 0 -q1;
%        -q2 q1 0];
% % q = [q1 q2 q3]';
% A = (q4^2 - dot(q,q))*eye(3) + 2*(q*q') - 2*q4*Q';
A = [q1^2 - q2^2 - q3^2 + q4^2, 2*(q1*q2 + q3*q4), 2*(q1*q3 - q2*q4);
2*(q1*q2 - q3*q4), -q1^2 + q2^2 - q3^2 + q4^2, 2*(q2*q3 + q1*q4);
2*(q1*q3 + q2*q4), 2*(q2*q3 - q1*q4), -q1^2 - q2^2 + q3^2 + q4^2];
% q = [q1 q2 q3 q4]';
% A = [q(1)^2-q(2)^2-q(3)^2+q(4)^2,...
%           2*(q(1)*q(2)+q(3)*q(4)),...
%           2*(q(1)*q(3)-q(2)*q(4));
%           2*(q(1)*q(2)-q(3)*q(4)),...
%          -q(1)^2+q(2)^2-q(3)^2+q(4)^2,...
%           2*(q(2)*q(3)+q(1)*q(4));
%           2*(q(1)*q(3)+q(2)*q(4)),...
%           2*(q(2)*q(3)-q(1)*q(4)),...
%          -q(1)^2-q(2)^2+q(3)^2+q(4)^2];
end