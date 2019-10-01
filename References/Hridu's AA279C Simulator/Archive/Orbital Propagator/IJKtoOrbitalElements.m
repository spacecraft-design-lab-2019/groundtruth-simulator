function oe = IJKtoOrbitalElements(state, GM)

oe = zeros(6, 1);

% Position and Velocity
rECI = state(1:3);
vECI = state(4:6);

% Specific Angular Momentum
h = cross(rECI, vECI);
W = h./norm(h);

% Inclination
i = atan2(norm(W(1:2)), W(3));

% Right Ascension of the Ascending Node
Omega = atan2(W(1),-W(2));
if Omega < 0
    Omega = Omega + 2*pi;
end

% Semi-Major Axis
a = (2/norm(rECI) - norm(vECI)^2/GM)^-1;

% Eccentricity
p = norm(h)^2/GM;
e = sqrt(1 - p/a);

% True Anomaly
n = sqrt(GM/a^3);
E = atan2((dot(rECI, vECI)/(n*a^2)), (1 - norm(rECI)/a));
if E < 0
    E = E + 2*pi;
end

nu = atan2(a*sqrt(1-e^2)*sin(E),a*(cos(E)-e));
if nu < 0
    nu = nu + 2*pi;
end

% Argument of Periapsis
u = atan2((rECI(3)/sin(i)), (rECI(1)*cos(Omega) + rECI(2)*sin(Omega)));
if u < 0
    u = u + 2*pi;
end

omega = u - nu;
if omega < 0
    omega = omega + 2*pi;
end

oe(1) = a;
oe(2) = e;
oe(3) = i;
oe(4) = Omega;
oe(5) = omega;
oe(6) = nu;

end