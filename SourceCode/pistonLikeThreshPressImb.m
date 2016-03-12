function [pistonLikeThreshPress,oilLayerExist] = pistonLikeThreshPressImb(geometry,half_angles,...
    rec_angle,adv_angle,inscribedR,G_factor,sig_ow)

global Pc_drain_max

hingeAngle = zeros(1,4); b = zeros(1,4); alfa = zeros(1,4);
if geometry == 1
    pistonLikeThreshPress = 2*sig_ow*cos(adv_angle) / inscribedR;
    oilLayerExist = nan;
else
    %% Piston like section
    if geometry == 2
        nc = 3;
    else
        nc = 4;
    end
    nominator = 0;
    for ii = 1:4
        if ~isnan(half_angles(ii))
            nominator = nominator + cos(rec_angle + half_angles(ii));
        end
    end
    maxAdvAngle = acos ((-4*G_factor*nominator)/...
        ((inscribedR*Pc_drain_max/sig_ow) - cos(rec_angle)+4*nc*G_factor*sin(rec_angle)));
    rpd = sig_ow/Pc_drain_max;
    if adv_angle <= maxAdvAngle
        rp1 = 2*rpd;
        rp2 = rpd;
        while abs(rp2 - rp1) > 10^-10
            rp1 = rp2;
            hingeAngle = acos((rpd / rp1)*cos(rec_angle + half_angles)) - half_angles;
            for ii = 1:4
                if ~isnan(half_angles(ii))
                    theta(ii) = min(hingeAngle(ii) , adv_angle);
                    E0(ii) = pi/2 - theta(ii) - half_angles(ii);
                    E1(ii) = cos(theta(ii) + half_angles(ii)) / sin(half_angles(ii));
                    E2(ii) = cos(theta(ii) + half_angles(ii))*cos(theta(ii)) / sin(half_angles(ii));
                end
            end
            rp2 = ((inscribedR^2 / 4 / G_factor) + rp1^2 * (sum(E0) - sum(E2))) /...
                (2*rp1 * sum(E0) + cos(adv_angle)*((inscribedR/2/G_factor) - 2*rp1*sum(E1)));
         end
        pistonLikeThreshPress = sig_ow / rp2;
        oilLayerExist = nan;
    elseif adv_angle > maxAdvAngle && adv_angle < pi/2 + max(half_angles)
        pistonLikeThreshPress = 2 * sig_ow * cos(adv_angle) / inscribedR;
        oilLayerExist = nan;
    elseif adv_angle >= pi/2 + max(half_angles)
        pistonLikeThreshPress = -throatThresholdPressPiston (pi - adv_angle,...
            geometry,G_factor,inscribedR,half_angles);
        oilLayerExist = 1;
    end
end