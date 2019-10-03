function M = GravTorques(mu,R,c,I)
cx = c(:,1);
cy = c(:,2);
cz = c(:,3);
Ix = I(1,1);
Iy = I(2,2);
Iz = I(3,3);
Mx = 3*mu./R.^3 * (Iz - Iy) .* cy.*cz;
My = 3*mu./R.^3 * (Ix - Iz) .* cz.*cx;
Mz = 3*mu./R.^3 * (Iy - Ix) .* cx.*cy;
M = [Mx, My, Mz];
end
