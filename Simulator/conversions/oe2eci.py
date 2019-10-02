# -*- coding: utf-8 -*-

def oe_2_eci(a, e, i, Omega, omega, nu, GM):
    """
    Function: oe_2_eci

    Converts orbital elements to ECI position and velocity state vector.

    Inputs:
        a:
        e: eccentricity [float64]
        i:
        Omega:
        omega:
        nu: polar angle (radians) [float64]
        GM:
        w: angular velocity (radians/s) np.array[3x1][float64]

    Outputs:
        q_dot: rate of change of quaternion np.array[4x1][float64]
    """
    #Solving for Eccentric Anomaly
    k = np.sqrt((1-e)/(1+e))
    E = 2*np.arctan(k*np.tan(nu/2))

    #Solving for Position in Perifocal Coordinates
    rPQW = np.array([a*(np.cos(E)-e), a*np.sqrt(1-e**2)*np.sin(E), 0])

    #Solving for Position in ECI Coordinates
    ROmega = np.array([[np.cos(-Omega), np.sin(-Omega), 0],
                        [-np.sin(-Omega), np.cos(-Omega), 0],
                        [0, 0, 1]])
    Ri = np.array([[1, 0, 0],
                    [0, np.cos(-i), np.sin(-i)],
                    [0, -np.sin(-i), np.cos(-i)]])
    Romega = np.array([[np.cos(-omega), np.sin(-omega), 0],
                        [-np.sin(-omega), np.cos(-omega), 0],
                        [0, 0, 1]])
    rotation = np.dot(np.dot(ROmega,Ri),Romega)
    rECI = np.dot(rotation,rPQW)

    #Solving for Velocity in Perifocal Coordinates
    n = np.sqrt(GM/a**3)
    factor = (a*n)/(1-e*np.cos(E))
    vPQW = np.dot(factor,np.array([-np.sin(E), np.sqrt(1-e**2)*np.cos(E), 0]))

    #Solving for Velocity in ECI Coordinates
    vECI = np.dot(rotation*vPQW)
    state = np.array([rECI, vECI])
    return state

"""
CONVERT TO PYTHON:


function state = calcPositionAndVelocity(a, e, i, Omega, omega, nu, GM)

% Solving for Eccentric Anomaly
k = sqrt((1-e)/(1+e));
E = 2*atand(k*tand(nu/2));

% Solving for Position in Perifocal Coordinates
rPQW = [a*(cosd(E) - e); a*sqrt(1 - e^2)*sind(E); 0];

% Solving for Position in ECI Coordinates
ROmega = [cosd(-Omega), sind(-Omega), 0;
    -sind(-Omega), cosd(-Omega), 0;
    0,            0,            1];
Ri =     [1,            0,            0;
    0,            cosd(-i),     sind(-i);
    0,           -sind(-i),     cosd(-i)];
Romega = [cosd(-omega), sind(-omega), 0;
    -sind(-omega), cosd(-omega), 0;
    0,            0,            1];
rotation = ROmega*Ri*Romega;
rECI = rotation*rPQW;

% Solving for Velocity in Perifocal Coordinates
n = sqrt(GM/a^3);
factor = (a*n)/(1 - e*cosd(E));
vPQW = factor.*[-sind(E); (sqrt(1 - e^2))*cosd(E); 0];

% Solving for Velocity in ECI Coordinates
vECI = rotation*vPQW;
state = [rECI; vECI];
end


"""
