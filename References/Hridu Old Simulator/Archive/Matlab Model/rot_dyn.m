function s_dot = rot_dyn(t, s, I, torque, type, w_r, w_r_dot, I_r, r)
n = length(s);
s_dot = zeros(n,1);

if strcmp(type, 'quaternions')
    s_dot(1:4) = calc_q_dot(s(1:4), s(5:7));
elseif strcmp(type, 'euler_angles')
    s_dot(1:3) = calc_e_dot(s(1:3), s(4:6));
end

s_dot(n-2:n) = calc_w_dot(s(n-2:n), I, torque, w_r, w_r_dot, I_r, r);
end