%how inaccurate can sun direction be and be accurate to 0.2 degrees?
true = [1;0;0];
n = 10000;
error = zeros(n,1);
noise = 0.002;
for i=1:n
    meas = true + noise*randn(size(true));
    meas = meas/norm(meas);
    error(i) = (180/pi)*acos(dot(true, meas));
end

figure, grid on, hold on
histogram(error)
xlabel('Error, degrees')