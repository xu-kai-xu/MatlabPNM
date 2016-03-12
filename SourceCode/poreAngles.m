function [minAngle , secondAngle , thirdAngle , forthAngle] = ...
    poreAngles (shapeFactor)
%%
% This function calculates half angle for each elemnt. circular elemnts have 
% four zero half angle, squares have four half angles equal to "pi/4" and 
% three half angle for triangular elemnts are calculated by Patzek method.

global n_p

for ii = 1:n_p
    %% Square Section
    if shapeFactor(ii) == 0.0625
        minAngle(ii) = pi / 4;
        secondAngle(ii) = pi / 4;
        thirdAngle(ii) = pi / 4;
        forthAngle(ii) = pi / 4;

    %% Circular Section
    elseif shapeFactor(ii) == 1/(4*pi);
        minAngle(ii) = nan;
        secondAngle(ii) = nan;
        thirdAngle(ii) = nan;
        forthAngle(ii) = nan;
    %% Triangular Section
    else
        bethaMin = atan(2*cos((acos(-12*sqrt(3)*shapeFactor(ii))/3)+(4*pi/3))/sqrt(3));
        bethaMax = atan(2*cos(acos(-12*sqrt(3)*shapeFactor(ii))/3)/sqrt(3));
        betha = bethaMin + rand*(bethaMax - bethaMin);

        alfa = -0.5 * betha + 0.5 * asin ((tan(betha) + 4*shapeFactor(ii)).*sin(betha)/...
            (tan(betha) - 4*shapeFactor(ii)));
        gama = pi /2 -(alfa + betha);
        
        sortedAngles = sort([alfa , betha , gama]);
        
        minAngle(ii) = sortedAngles(1) ;
        secondAngle(ii) = sortedAngles(2);
        thirdAngle(ii) = sortedAngles(3);
        forthAngle(ii) = nan;
    end
end