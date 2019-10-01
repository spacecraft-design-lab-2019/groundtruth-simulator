function state = calcPositionAndVelocity(a, e, i, Omega, omega, nu, GM)

    % Solving for Eccentric Anomaly
    k = sqrt((1-e)/(1+e));
    E = 2*atand(k*tand(nu/2));
    
    % Solving for Position in Perifocal Coordinates
    rPQW = [a*(cosd(E) - e); a*sqrt(1 - e^2)*sind(E); 0];
    
    % Solving for Position in ECI Coordinates
    ROmega = [cosd(-Omega), sind(-Omega), 0;
        -sind(-Omega), cosd(-Omega), 0;
        0,            0,            1];
    Ri =     [1,            0,            0;
        0,            cosd(-i),     sind(-i);
        0,           -sind(-i),     cosd(-i)];
    Romega = [cosd(-omega), sind(-omega), 0;
        -sind(-omega), cosd(-omega), 0;
        0,            0,            1];
    rotation = ROmega*Ri*Romega;
    rECI = rotation*rPQW;
    
    % Solving for Velocity in Perifocal Coordinates
    n = sqrt(GM/a^3);
    factor = (a*n)/(1 - e*cosd(E));
    vPQW = factor.*[-sind(E); (sqrt(1 - e^2))*cosd(E); 0];
    
    % Solving for Velocity in ECI Coordinates
    vECI = rotation*vPQW;
    state = [rECI; vECI];
end