function [phi,theta,psi] = quat2EA313(q)
q1t = q(:,1);
q2t = q(:,2);
q3t = q(:,3);
q4t = q(:,4);

phi = unwrap(atan2(2*q1t.*q3t - 2*q4t.*q2t,2*q2t.*q3t+ 2*q4t.*q1t));
theta = acos(q3t.^2 - q2t.^2 - q1t.^2 + q4t.^2);
psi = unwrap(atan2(2*q1t.*q3t+2*q4t.*q2t,-2*q2t.*q3t+ 2*q4t.*q1t));

end