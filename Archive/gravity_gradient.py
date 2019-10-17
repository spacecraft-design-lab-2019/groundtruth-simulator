def gravityGradientTorque(r_sat, R_eci2principal, I, GM):
	    """
    Function: gravityGradientTorque

    Inputs:
        r_sat: position vector of satellite in ECI
        R_eci2principal: rotation matrix from eci to principal axes
        I: moment of inertia matrix (in principal axes) (diagonal)
        GM: gravitational parameter of attracting body
    Outputs:
        M: moment due to gravity gradient
    """

    
    # R = np.linalg.norm(r_sat)
    # c = R_eci2principal @ r_sat
    
    # Ix = I[0,0]
    # Iy = I[1,1]
    # Iz = I[2,2]
    # cx = c[0]
    # cy = c[1]
    # cz = c[2]
    
    # M = np.zeros((3,))
    # M[0] = 3*GM/R**3 * (Iz - Iy)*cy*cx
    # M[1] = 3*GM/R**3 * (Ix - Iz)*cz*cx
    # M[2] = 3*GM/R**3 * (Iy - Ix)*cx*cy
    
    # return M

