function [ R ] = calc_R( e )
%Calculate Rotation matrix, given Euler angles
%INPUT: euler angles, roll pitch yaw
%OUTPUT: rotation matrix
RX = [1,0,0; 
      0,cos(e(1)), -1*sin(e(1));
      0,sin(e(1)),    cos(e(1))];
RY = [cos(e(2)), 0, sin(e(2));
      0,         1, 0;
   -1*sin(e(2)), 0, cos(e(2))];
RZ = [cos(e(3)), -1*sin(e(3)), 0;
      sin(e(3)),    cos(e(3)), 0;
      0, 0,                    1];

R = RZ*RY*RX;
end

