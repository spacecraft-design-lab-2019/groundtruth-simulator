#def calc_e_dot(e,w):
#    """
#    Function: calc_e_dot
#
#    Calculates derivative of euler angles.
#    Input:
#    e: euler angle vector (radians) np.array[3x1][float64]
#    w: angular velocity vector (radians/s) np.array[3x1][float64]
#    Output:
#    e_dot: rate of change of eular angles np.array[1x3][float64]
#    """
#    #propagate actual dynamics
#    tan_the = np.tan(e[1])
#    if np.abs(tan_the) > 300: #make sure tan(theta) is well defined
#        tan_the = 300*np.sign(tan_the)
#    cos_the = np.cos(e[1])
#    if np.abs(cos_the) < 10e-4: #make sure cos(theta) is well defined
#        cos_the = 10e-4*np.sign(cos_the)
#    A = np.array([[1, tan_the*np.sin(e[0]), tan_the*np.cos(e[0])],
#                [0, np.cos(e[0]), -1*np.sin(e[0])],
#                [0, np.sin(e[0])/cos_the, np.cos(e[0])/cos_the]])
#    return np.dot(A,w)