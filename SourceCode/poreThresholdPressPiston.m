function poreThresholdPressure = poreThresholdPressPiston (recedingAngle,...
    poreGeometry,shapeFactor,poreRadius,halfAngles)
% Calculating element thrshold pressure according to Patzek's approach

global sig_ow 
if poreGeometry == 1 %circular
    poreThresholdPressure = 2*sig_ow*cos(recedingAngle)/poreRadius;
else
    nominator = 0;
    for jj = 1:4
        if ~isnan(halfAngles(jj)) && halfAngles(jj) < pi/2 - recedingAngle
            E2 = cos(recedingAngle + halfAngles(jj)) *...
                cos(recedingAngle) / sin(halfAngles(jj));
            E0 = pi / 2 - recedingAngle - halfAngles(jj);
            nominator = nominator +  (E2 - E0);
        end
    end
    poreThresholdPressure = (sig_ow / poreRadius)*cos(recedingAngle)*...
        (1+sqrt(1 - (4*shapeFactor*nominator)/(cos(recedingAngle)^2)));
end