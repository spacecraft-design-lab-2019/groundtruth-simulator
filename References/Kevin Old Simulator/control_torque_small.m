function M_c = control_torque_small(I, w, Mp, A_E)

f = 1000;

kp = f^2./diag(I);
kd = 2*sqrt(diag(I).*(Mp + kp));
M_c(1) = -kd(1)*w(1) - kp(1)*(A_E(2,3));
M_c(2) = -kd(2)*w(2) - kp(2)*(A_E(3,1));
M_c(3) = -kd(3)*w(3) - kp(3)*(A_E(1,2));
M_c = saturated(M_c, -500, 500);
end

function Msat = saturated(M, lo, hi)
Msat = min(max(M, lo),hi);
end