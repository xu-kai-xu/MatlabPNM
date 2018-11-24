classdef Element
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % properties set from the link and node input files
        index
        radius
        shapeFactor
        volume
        clayVolume
        isInlet
        isOutlet
        
        % Calculated properties
        conductance %Element conductance
        geometry % geometrical shape of the element
        halfAngle1
        halfAngle2
        halfAngle3
        halfAngle4
        area
        waterCrossSectionArea
        oilCrossSectionArea
        gasCrossSectionArea 
        oilArea
        waterArea
        
        %It's better to define the phase conductances and the area of each
        %fluid as a structure for this each element
        waterConductance
        oilConductance 
        gasConductance
        fluidConductances
        fluidCrossSectionAreas

        receedingContactAngle=20;
        advandingContactAngle
        
        waterPressure
        oilPressure
        gasPressure
        thresholdPressure
        
     
    end
    
    methods
        function ThresholdPressure = calculateThresholdPressurePistonLike (obj, sig_ow)
         % calculateThresholdPressurePistonLike Summary of this method goes here
         % Detailed explanation goes here         
             if strcmp(obj.geometry , 'Circle')== 1
                 ThresholdPressure = 2*sig_ow*cos(obj.receedingContactAngle)/obj.radius;
             else
                 nominator = 0;
                 halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                 for jj = 1:4
                     if ~isnan(halfAngles(jj)) && halfAngles(jj) < pi/2 - obj.receedingContactAngle
                         E2 = cos(obj.receedingContactAngle + halfAngles(jj)) *...
                             cos(obj.receedingContactAngle) / sin(halfAngles(jj));
                         E0 = pi / 2 - obj.receedingContactAngle - halfAngles(jj);
                         nominator = nominator +  (E2 - E0);
                     end
                 end
                 ThresholdPressure = (sig_ow / obj.radius)*...
                     cos(obj.receedingContactAngle)*(1+sqrt(1 -(4*obj.shapeFactor*...
                     nominator)/(cos(obj.receedingContactAngle)^2)));
             end
         end
        
         %% Conductance Calculation
         function [waterArea, waterConductance] = calculateWaterConductance(obj, network)           
                curvatureRadius = network.sig_ow / obj.thresholdPressure;
                b = zeros(1,4); % Meniscus-Apex distance
                if strcmp(obj.geometry , 'Circle')== 1
                    waterArea = 0;
                    waterConductance = 0;
                else
                    dimenlessAreaCorner = zeros(1,4);
                    cornerShapeFactor = zeros(1,4);
                    scaledDimenlessConducCorner = zeros(1,4);
                    dimenlessConducCorner = zeros(1,4);
                    halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                    for jj = 1:4
                        if ~isnan(halfAngles(jj))
                            b(jj) = (curvatureRadius * cos(obj.receedingContactAngle +...
                                halfAngles(jj)) / sin(halfAngles(jj)));
                            
                            if (halfAngles(jj) + obj.receedingContactAngle) == pi/2
                                dimenlessAreaCorner(jj) = sin(halfAngles(jj))*cos(halfAngles(jj));
                                cornerShapeFactor(jj) = dimenlessAreaCorner(jj) /...
                                    (4 * (1 + sin(halfAngles(jj)))^2);
                            else
                                dimenlessAreaCorner(jj) = (sin(halfAngles(jj)) /...
                                    cos(obj.receedingContactAngle + halfAngles(jj)))^2 *...
                                    ((cos(obj.receedingContactAngle)*...
                                    cos(obj.receedingContactAngle + halfAngles(jj)) /...
                                    sin(halfAngles(jj)))+ obj.receedingContactAngle +...
                                    halfAngles(jj) - pi/2);
                                cornerShapeFactor(jj) = dimenlessAreaCorner(jj) /...
                                    (4 * (1 - (sin(halfAngles(jj)) / cos(halfAngles(jj)+...
                                    obj.receedingContactAngle))*...
                                    (obj.receedingContactAngle + halfAngles(jj) - pi/2))^2);
                            end
                            
                            scaledDimenlessConducCorner(jj) = -15.1794 * cornerShapeFactor(jj)^2 +...
                                7.6307 * cornerShapeFactor(jj) - 0.53488;
                            
                            dimenlessConducCorner(jj) = exp((scaledDimenlessConducCorner(jj) +...
                                0.02 * sin(halfAngles(jj) - pi/6)) / (1/4/pi - cornerShapeFactor(jj))^(7/8)...
                                / cos(halfAngles(jj) - pi/6)^0.5) * dimenlessAreaCorner(jj)^2;
                        end
                    end
                    waterArea = sum(b.^2 .* dimenlessAreaCorner)
                    waterConductance = sum(2 * b.^4 .* dimenlessConducCorner / network.waterViscosity);
                end
         end
         function [oilArea, oilConductance] = calculateOilConductance(obj, network)    
             oilArea = obj.area - obj.waterArea;
             if strcmp(obj.geometry , 'Circle')== 1
                 oilConductance = (oilArea / obj.area)*0.5 * obj.area^2 * obj.shapeFactor /network.oilViscosity;
             elseif strcmp(obj.geometry , 'Triangle')== 1
                 oilConductance = (oilArea) *3*obj.area^2* obj.shapeFactor/network.oilViscosity/5;
             elseif strcmp(obj.geometry , 'Square')== 1
                 oilConductance = (oilArea / obj.area)*0.5623 * obj.area^2 * obj.shapeFactor /network.oilViscosity;
             end              
         end
    end

end

