# -*- coding: utf-8 -*-
"""
Function: accelPointMass

Calculates the acceleration due to a point mass. Returns the derivative of the
current state vector.
"""

def accelPointMass(state, GM):
    
    
"""
CONVERT TO PYTHON

function statedot = f(state, GM)
    statedot = zeros(6,1);
    statedot(1:3) = state(4:6);
    statedot(4:6) = -GM*state(1:3)./norm(state(1:3))^3;
end

"""