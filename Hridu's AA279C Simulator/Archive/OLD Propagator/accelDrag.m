function f_drag = accelDrag(state)

X = load('constants.mat');

r = state(1:3);
v = state(4:6);

B = X.Cd_drag*X.A_drag/X.M_drag;
h = norm(r);
p = X.p0_drag*exp(-(h - X.h0_drag)/X.H_drag);

vrel = v - cross([0; 0; X.omega_e], r);

f_drag = -0.5.*B.*p.*(norm(vrel)^2).*vrel./norm(vrel);
end