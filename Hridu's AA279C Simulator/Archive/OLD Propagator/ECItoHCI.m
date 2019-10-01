function rHCI  = ECItoHCI(obq, rECI)
rotation = [1,   0,          0;
            0,   cosd(-obq),  sind(-obq);
            0,  -sind(-obq), cosd(-obq)];
rHCI = rotation*rECI;
end