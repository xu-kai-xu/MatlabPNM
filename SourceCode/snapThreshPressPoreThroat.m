function snapThreshPress = snapThreshPressPoreThroat(geometry,half_angles,...
    rec_angle,adv_angle,inscribedR,sig_ow)

global Pc_drain_max

% hingeAngle = zeros(1,4);
if geometry == 1
    snapThreshPress = nan;
else
    %% Snap off Section
    maxAdvAngle = pi/2 - min(half_angles);
    if geometry == 2
        nc = 3;
    elseif geometry == 3;
        nc = 4;
    end
    if adv_angle < maxAdvAngle
        rso = zeros(1,nc);
        for edge = 1:nc
            rso1 = 1; rso2 = 2;
            while abs(rso1 - rso2) > 10^-8
                rso1 = rso2;
                ii = edge;
                hingeAngle_ii = acos((sig_ow / rso1)*...
                    cos(rec_angle+half_angles(ii))/Pc_drain_max) - half_angles(ii);
                theta_i = min(hingeAngle_ii , adv_angle);
                E1_i = cos(theta_i + half_angles(ii))/sin(half_angles(ii));
                if ii == nc
                    jj = 1;
                    hingeAngle_jj = acos((sig_ow / rso1)*...
                        cos(rec_angle+half_angles(jj))/Pc_drain_max) - half_angles(jj);
                    theta_j = min(hingeAngle_jj , adv_angle);
                    E1_j = cos(theta_j + half_angles(jj))/sin(half_angles(jj));
                    rso2 = inscribedR *(cot(half_angles(ii))+cot(half_angles(jj)))/(E1_i+E1_j);
                    rso(edge) = rso2;
                else
                    jj = edge + 1;
                    hingeAngle_jj = acos((sig_ow / rso1)*...
                        cos(rec_angle + half_angles(jj)) / Pc_drain_max) - half_angles(jj);
                    theta_j = min(hingeAngle_jj , adv_angle);
                    E1_j = cos(theta_j + half_angles(jj))/sin(half_angles(jj));
                    rso2 = inscribedR *(cot(half_angles(ii))+cot(half_angles(jj)))/(E1_i+E1_j);
                    rso(edge) = rso2;
                end
            end
        end   
        snapThreshPress = sig_ow / min(rso);   
    elseif adv_angle == maxAdvAngle
        snapThreshPress = 0;
    elseif adv_angle > maxAdvAngle && adv_angle < pi - min(half_angles)
        snapThreshPress = Pc_drain_max*cos(adv_angle + min(half_angles))/...
            cos(rec_angle + min(half_angles));
    elseif adv_angle >= pi - min(half_angles)
        snapThreshPress = -Pc_drain_max/cos(rec_angle + min(half_angles));
    end 
end