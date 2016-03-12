function [oil_conductance, water_conductance, A_oil, A_water] =...
    conduct_imb(Pc_threshold,rec_angle,adv_angle,half_angles,geometry,...
    total_A, layerExist , single_phase_conduct,G_factor, R_ins , PcMax , PcLayerCorner)

% Piri's Approach
% PcMax is the Pc which the layer formed in.

global water_viscosity oil_viscosity sig_ow Pc_drain_max
if geometry == 1
    A_oil = 0;
    A_water = total_A;
    oil_conductance = 0;
    water_conductance = single_phase_conduct;
else
    if geometry == 2
        nc = 3;
    elseif geometry == 3
        nc = 4;
    end
    
    if isnan(layerExist) % If no oil layer exist
    A_oil = 0;
    A_water = total_A;
    oil_conductance = 0;
    water_conductance = single_phase_conduct;
    else % If oil layer exist
        R = sig_ow / abs(Pc_threshold);
        A_outer = zeros(1,4);
        A_iner = zeros(1,4);
        g_iner = zeros(1,4);
        A_layer = zeros(1,4);
        g_layer = zeros(1,4);
        for ii = 1:nc
            if ~isnan(PcLayerCorner(ii)) % if the layer exist in the corner
                Rmin = sig_ow / Pc_drain_max;
                bi = Rmin * (cos(rec_angle + half_angles(ii)) / sin(half_angles(ii)));
                thetaHingeRec = acos(bi * sin(half_angles(ii)) / R) - half_angles(ii);
                if thetaHingeRec > adv_angle
                    thetaHingeRec = adv_angle;
                end
                % Area of corner water
                if half_angles(ii) + thetaHingeRec == pi/2
                    A_iner(ii) = (R*cos(thetaHingeRec + half_angles(ii))/...
                        sin(half_angles(ii)))^2 * sin(half_angles(ii))*cos(half_angles(ii));
                else
                    A_iner(ii) = R^2*(cos(thetaHingeRec)*(cot(half_angles(ii))*...
                        cos(thetaHingeRec)-sin(thetaHingeRec))+thetaHingeRec + half_angles(ii) - pi/2);
                end
                % Conductance of corner water
                if thetaHingeRec <= pi/2 - half_angles(ii)
                    phi1 = pi/2 - half_angles(ii) - thetaHingeRec;
                    phi2 = cot(half_angles(ii))*cos(thetaHingeRec) - sin(thetaHingeRec);
                    phi3 = (pi/2 - half_angles(ii))*tan(half_angles(ii));
                    g_iner(ii) = (A_iner(ii)^2*...
                        (1-sin(half_angles(ii)))^2*(phi2*cos(thetaHingeRec)-phi1)*phi3^2)/...
                        (12*water_viscosity*sin(half_angles(ii))^2*(1-phi3)^2*(phi2 + phi1)^2);
                elseif thetaHingeRec > pi/2 - half_angles(ii)
                    g_iner(ii) = (A_iner(ii)^2*tan(half_angles(ii))*...
                        (1-sin(half_angles(ii)))^2*phi3^2)/(12*water_viscosity*...
                        sin(half_angles(ii))^2*(1-phi3)*(1+phi3)^2);
                end
                
                Rmin_outer = sig_ow / abs(PcMax);
%                 bo = Rmin_outer * (cos(pi - adv_angle+half_angles(ii)) / sin(half_angles(ii)));
                bo = sig_ow * (cos(pi - adv_angle+half_angles(ii)) / sin(half_angles(ii))) / abs(Pc_threshold);
                thetaHingeAdv = pi - adv_angle;
%                 thetaHingeAdv = pi - (acos(bo*sin(half_angles(ii)) / R) - half_angles(ii));
                % Area of corner fluid(oil+water)
                if half_angles(ii) + thetaHingeAdv == pi/2
                    A_outer(ii) = (R*cos(thetaHingeAdv + half_angles(ii))/...
                        sin(half_angles(ii)))^2*sin(half_angles(ii))*cos(half_angles(ii));
                else
                    A_outer(ii) = R^2*(cos(thetaHingeAdv)*(cot(half_angles(ii))*...
                        cos(thetaHingeAdv)-sin(thetaHingeAdv))+thetaHingeAdv+half_angles(ii)-pi/2);
                end                
                % Area of oil layer
                A_layer(ii) = A_outer(ii) - A_iner(ii);
                % Conductance of oil layer
                phi3_layer = (pi/2 - half_angles(ii))*tan(half_angles(ii));
                g_layer(ii) = (A_layer(ii)^3*...
                    (1-sin(half_angles(ii)))^2 * tan(half_angles(ii)) * phi3_layer^2)/...
                    (12*oil_viscosity*A_outer(ii)*sin(half_angles(ii))^2*...
                    (1-phi3_layer)*(1+phi3_layer-(1-phi3_layer)*...
                    sqrt(A_iner(ii)/A_outer(ii))));
            end
        end
        
        % Center water area and conductance
        A_water_center = total_A - sum(A_outer);
        if geometry == 1
            g_water_center = 0.5*G_factor*A_water_center^2/water_viscosity;
        elseif geometry == 2
            g_water_center = 3*R_ins^2*A_water_center/20/water_viscosity;
        elseif geometry == 3
            g_water_center = 0.5623*G_factor*A_water_center^2 /water_viscosity;
        end
        A_water = total_A - sum(A_layer);
        water_conductance = g_water_center + sum(g_iner);
        
        A_oil = sum(A_layer);
        oil_conductance = sum(g_layer);
    end
end