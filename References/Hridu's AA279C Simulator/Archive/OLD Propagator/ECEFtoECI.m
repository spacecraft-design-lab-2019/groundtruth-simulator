function rECI  = ECEFtoECI(GMST, rECEF)
rotation = [cosd(-GMST), sind(-GMST), 0;
           -sind(-GMST), cosd(-GMST), 0;
            0,          0,          1];
rECI = rotation*rECEF;
end