function e = quat2eul(q_in)
%translate quaternion to euler angles
%v1.0, Apr 17 2019
%assumes zyx, as present in the existing kinematics model
q1 = q_in(1); q2 = q_in(2); q3 = q_in(3); q4 = q_in(4);
a =  2*(q2*q3 + q1*q4);
b = q1*q1 + q2*q2 - q3*q3 - q4*q4;
c = -2*(q2*q4 - q1*q3);
d =  2*(q3*q4 + q1*q2);
e = q1*q1 - q2*q2 - q3*q3 + q4*q4;
%these trig functions will change if we change rotation order
rol = atan2(a, b);
pit = asin(c);
yaw = atan2(d,e);

e = [rol, pit, yaw];
end