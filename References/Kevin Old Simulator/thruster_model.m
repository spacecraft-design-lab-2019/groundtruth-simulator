function T_real = thruster_model(T_des, dt, prev_T_real)

%quantize input
ntBP = numerictype(1,8,4);
y = quantize(fi(T_des),ntBP);

%limit rate
R = 6; %Rise Slew Rate Max
F = -6; %Fall Slew Rate Max
y = rate_limiter(y, dt, prev_T_real, R, F);

%limit_thrust
maxt = 4.5;
y = thrust_limiter(y, maxt);

% add noise
y = add_noise(y,dt); 

mint = 0.5;
y = dead_zone(y,mint);

T_real = y;
end

function y = dead_zone(y, mint)
if y < mint
    y = 0;
end
end

function y = add_noise(y,dt)
mu = 0;
sig = 0.0006/sqrt(dt);
y = y + mvnrnd(mu,sig);
end

function y = thrust_limiter(y, maxt)
y = min(y, maxt);

end

function y = rate_limiter(T_des, dt, prev_T_real, R, F)
u = T_des;
y = prev_T_real;

rate = (u - y)/dt;

if rate > R
    y = dt*R + y;
elseif rate < F
    y = dt*F + y;
else
    y = u;    
end

end