function w_real = wheel_model(Tor_des, dt, prev_Tor)
%quantize input
ntBP = numerictype(1,8,4);
y = quantize(fi(Tor_des),ntBP);

% add noise
y = add_noise(y,dt); 

%saturation
minT = -3;
maxT = 3;
y = saturation(y, minT, maxT);

w_real = y;
end

function y = saturation(y, minT, maxT)
y = max(min(y, maxT), minT);
end

function y = add_noise(y,dt)
mu = 0;
sig = 0.002/sqrt(dt);
y = y + mvnrnd(mu,sig);
end