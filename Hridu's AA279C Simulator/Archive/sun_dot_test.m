%test direction and incident angle
sat2sun_eci = zeros(3,1);
n = 20;
x = linspace(-1, 1, n);
y = linspace(-1, 1, n);
z = linspace(-1, 1, n);
test = [1,-1,0, 0,0, 0;
        0, 0,1,-1,0, 0;
        0, 0,0, 0,1,-1];
    
sun_meas = zeros(6,1);
for a=1:n
    for b=1:n
        for c = 1:n
            sun = [x(a); y(b); z(c)];
            sun = sun/norm(sun);
            
            for i=1:6
                sun_meas(i) = dot(test(:,i), sun);
                if(sun_meas(i) < 0)
                    sun_meas(i) = 0;
                end
            end
            sun_meas
        end
    end
end