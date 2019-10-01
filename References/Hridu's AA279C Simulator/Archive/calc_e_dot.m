function [ e_dot ] = calc_e_dot( e, w )
%INPUT: euler angles, angular velocity
%OUTPUT: rate of change of euler angles
%propagate actual dynamics
    tan_the = tan(e(2));
    if(abs(tan_the) > 300) %make sure tan(theta) is well defined
        tan_the = 300*sign(tan_the);
        disp('tan x')
    end
    cos_the = cos(e(2));
    if(abs(cos_the) < 10^-4) %make sure cos(theta) is well defined
        cos_the = 10^-4*sign(cos_the);
        disp('cos x')
    end
    A = [1, tan_the*sin(e(1)), tan_the*cos(e(1));
         0, cos(e(1)),              -1*sin(e(1));
         0, sin(e(1))/cos_the,  cos(e(1))/cos_the];
    e_dot = A*w;
end

