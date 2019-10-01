# -*- coding: utf-8 -*-
"""
Function: calc_w_dot

Uses Euler equations to calculate derivative of angular velocity vector.
"""

def calc_w_dot():
    
    

"""
CONVERT TO PYTHON

function w_dot = fcn(w, tor, I)

w_dot = -I\(Skew(w)*I*w-tor);

end


function S = Skew(v)
S = [0 -v(3) v(2)
    v(3) 0 -v(1)
    -v(2) v(1) 0];
end

"""