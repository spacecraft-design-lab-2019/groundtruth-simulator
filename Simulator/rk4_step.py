# -*- coding: utf-8 -*-
def rk4_step(f, t0, x0, h):
    """
    Function: rk4_step
        Uses rk4 to integrate a single step. (vectorized)
        
    Inputs:
        f:  derivative of state vector
        t0: initial time
        x0: initial state vector
        h:  step size (in seconds)
    Outputs:
        t1: updated time
        x1: updated state vector
    """
    
    k1 = h * f(t0,      x0)
    k2 = h * f(t0+h/2,  x0+k1/2)
    k3 = h * f(t0+h/2,  x0+k2/2)
    k4 = h * f(t0+h,    x0+k3)
    x1 = x0 + (k1 + 2*k2 + 2*k3 + k4)/6
    t1 = t0 + h
    return t1, x1