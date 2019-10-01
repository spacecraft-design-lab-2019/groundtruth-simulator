function rECI  = HCItoECI(obq, rHCI)
rotation = [1,   0,          0;
            0,   cosd(obq), sind(obq);
            0,  -sind(obq), cosd(obq)];
rECI = rotation*rHCI;
end