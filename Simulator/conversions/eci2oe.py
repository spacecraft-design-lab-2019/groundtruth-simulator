# -*- coding: utf-8 -*-
"""
Function: eci2oe

Converts ECI position and velocity to orbital elements.
"""

def eci2oe(state, GM):
    
    
"""
CONVERT TO PYTHON:
    
function oe = eci2oe(state, GM)

r_ijk = state(1:3);
v_ijk = state(4:6);
r = norm(r_ijk);
v = norm(v_ijk);

% Create all necessary vectors
hVec = cross(r_ijk, v_ijk);
h = norm(hVec);

nVec = cross([0, 0, 1], hVec);
n = norm(nVec);

eVec = (1/GM)*((v^2 - GM/r)*r_ijk - dot(r_ijk, v_ijk)*v_ijk);
e = norm(eVec);

% Compute the size of the orbit
mechEnergy = 0.5*v^2 - GM/r;

if e ~=1
    a = -GM/(2*mechEnergy);
    p = a*(1 - e^2);
else
    a = inf;
    p = h^2/GM;
end

% Compute the orientation of the orbit
i = acosd(hVec(3)/h);
Omega = acosd(nVec(1)/n);
omega = acosd(dot(nVec, eVec)/(n*e));
nu = acosd(dot(eVec, r_ijk)/(e*r));

% Place angles in the correct domains
if nVec(2) < 0
    Omega = 360 - Omega;
end
if eVec(3) < 0
    omega = 360 - omega;
end
if dot(r_ijk, v_ijk) < 0
    nu = 360 - nu;
end

% Account for any special cases
if i == 0 && e ~= 0  % Elliptical equatorial
    % Provide the longitude of periapsis (PI = Om + w)
    ang = acosd(eVec(1)/e);
    
    if eVec(2) < 0
        ang = 360 - ang;
    end
elseif i~= 0 && e == 0 % Circular inclined
    % Provide the argument of latitude (u = w + anom)
    ang = acosd(dot(nVec, r_ijk)/(n*r));
    
    if r_ijk(3) < 0
        ang = 360 - ang;
    end
elseif i == 0 && e == 0 % Circular equatorial
    % Provide the true latitude (lambda = Om + w + anom)
    ang = acosd(r_ijk(1)/r);
    
    if r_ijk(2) < 0
        ang = 360 - ang;
    end
else
    % Default output for ang
    ang = NaN;
end

oe(1) = a;
oe(2) = e;
oe(3) = i;
oe(4) = Omega;
oe(5) = omega;
oe(6) = nu;
oe(7) = ang;

oe(3:7) = deg2rad(oe(3:7));
end


"""