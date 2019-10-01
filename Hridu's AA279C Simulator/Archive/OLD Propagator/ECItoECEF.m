function rECEF  = ECItoECEF(GMST, rECI)
rotation = [cosd(GMST), sind(GMST), 0;
           -sind(GMST), cosd(GMST), 0;
            0,          0,          1];
rECEF = rotation*rECI;
end