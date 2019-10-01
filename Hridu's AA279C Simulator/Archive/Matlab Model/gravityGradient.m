function M = gravityGradient(I, GM, rECI, R_eci2principal)
% Expressed in principal axes. The c vector is the position vector of the
% satellite expressed in principal axes.

R = norm(rECI);
c = R_eci2principal*rECI;

Ix = I(1,1);
Iy = I(2,2);
Iz = I(3,3);
cx = c(1);
cy = c(2);
cz = c(3);

M = zeros(3,1);
M(1) = 3*GM/R^3 * (Iz - Iy)*cy*cz;
M(2) = 3*GM/R^3 * (Ix - Iz)*cz*cx;
M(3) = 3*GM/R^3 * (Iy - Ix)*cx*cy;
end