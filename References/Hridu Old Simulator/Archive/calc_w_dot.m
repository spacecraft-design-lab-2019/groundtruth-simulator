function w_dot = calc_w_dot(w, I, tor, w_r, w_r_dot, I_r, r)
    %w_dot without rotor
    w_dot(1, 1) = (tor(1) - (I(3,3)-I(2,2))*w(2)*w(3))/I(1,1);
    w_dot(2, 1) = (tor(2) - (I(1,1)-I(3,3))*w(3)*w(1))/I(2,2);
    w_dot(3, 1) = (tor(3) - (I(2,2)-I(1,1))*w(1)*w(2))/I(3,3);
    
    %w_dot with rotor, assume rotor speed is constant or w_r_dot is known
    wrx = I_r*w_r_dot*r(1) + I_r*w_r*(w(2)*r(3) - w(3)*r(2));
    wry = I_r*w_r_dot*r(2) + I_r*w_r*(w(3)*r(1) - w(1)*r(3));
    wrz = I_r*w_r_dot*r(3) + I_r*w_r*(w(1)*r(2) - w(2)*r(1));
    w_dot(1, 1) = (tor(1) - (I(3,3)-I(2,2))*w(2)*w(3) - wrx) /I(1,1);
    w_dot(2, 1) = (tor(2) - (I(1,1)-I(3,3))*w(3)*w(1) - wry) /I(2,2);
    w_dot(3, 1) = (tor(3) - (I(2,2)-I(1,1))*w(1)*w(2) - wrz) /I(3,3);
end