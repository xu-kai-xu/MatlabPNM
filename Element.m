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
        water_conductance
        oil_conductance         
        geometry % geometrical shape of the element
        halfAngle1
        halfAngle2
        halfAngle3
        halfAngle4
        area
<<<<<<< HEAD
<<<<<<< HEAD
        waterCrossSectionArea
        oilCrossSectionArea
        gasCrossSectionArea 

||||||| merged common ancestors
        oil_area
        water_area
=======
        waterCrossSectionArea
        oilCrossSectionArea
        gasCrossSectionArea 
        oilArea
        waterArea
>>>>>>> 61a63595e6d7105ef83e8af42fce616fbb52602f
||||||| ad1e304... Change conductance formula
        waterCrossSectionArea
        oilCrossSectionArea
        gasCrossSectionArea 

=======
        oil_area
        water_area
>>>>>>> parent of ad1e304... Change conductance formula
        
        %It's better to define the phase conductances and the area of each
        %fluid as a structure for this each element
        fluidConductances
        fluidCrossSectionAreas
        
        receedingContactAngle=20;
        advandingContactAngle
        
        waterPressure
        oilPressure
        gasPressure
        
     
    end
    
    methods
<<<<<<< HEAD
<<<<<<< HEAD
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
         function [waterCrossSectionArea, waterConductance] = calculateWaterConductance(obj, network) 
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
             cornerArea = zeros(1,4);
             cornerConductance = zeros(1,4);
             %Raduis of Curvature
             Rc = network.waterViscosity / obj.thresholdPressure;
             for jj = 1:4
                 if ~isnan(halfAngles(jj))
                     % Based on Piri_2005: eq A4 & A5 
                     if (halfAngles(jj) + obj.receedingContactAngle) == pi/2
                          cornerArea(jj) = (Rc * cos(obj.receedingContactAngle + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                              sin(halfAngles(jj)) * cos(halfAngles(jj));
                      else
                          cornerArea(jj) = Rc ^2 * (cos(obj.receedingContactAngle)*...
                              (cot(halfAngles(jj)) * cos(obj.receedingContactAngle) - sin(obj.receedingContactAngle))+ ...
                              obj.receedingContactAngle + halfAngles(jj) - pi/2);
                     end   
                     % Based on Piri_2005: eq B(10 - 15)
                     f = 1;
                     F1 = pi/2 - halfAngles(jj) - obj.receedingContactAngle;
                     F2 = cot(halfAngles(jj)) * cos(obj.receedingContactAngle) - sin(obj.receedingContactAngle); 
                     F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));
                     
                     if (obj.receedingContactAngle <= pi/2 - halfAngles(jj))                         
                         cornerConductance(jj) = (cornerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                             (F2 * cos(obj.receedingContactAngle) - F1) * F3 ^ 2) / ...
                             (12 * network.waterViscosity * ((sin(halfAngles(jj))) * ...
                             (1 - F3) * (F2 + f * F1))^ 2);
                     elseif (obj.receedingContactAngle > pi/2 - halfAngles(jj))
                         cornerConductance(jj) = (cornerArea(jj)^2 * tan(halfAngles(jj))* ...
                             (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                             (12 * network.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
%                      else
%                           cornerConductance(jj) = (A1^3 * (1 - sin(halfAngles(jj)))^2 * ...
%                              tan(halfAngles(jj)) * F3 ^ 2) / ...
%                               (12 * network.waterViscosity * (sin(halfAngles(jj)))^2 * ...
%                               (1 - F3) * (1 + f1*F3 - (1- f2*F3) * sqrt(A2/cornerArea(jj)))^ 2);
                     end
                 end
                 waterCrossSectionArea = sum(cornerArea);
                 waterConductance = sum(cornerConductance);
             end
         end
           
         function [oilCrossSectionArea, oilConductance] = calculateOilConductance(obj, network)    
             oilCrossSectionArea = obj.area - obj.waterCrossSectionArea;
             if strcmp(obj.geometry , 'Circle')== 1
                 oilConductance = (oilCrossSectionArea / obj.area)*0.5 * obj.area^2 * ...
                     obj.shapeFactor /network.oilViscosity;
             elseif strcmp(obj.geometry , 'Triangle')== 1
                 oilConductance = oilCrossSectionArea * 3 * obj.radius^2 / network.oilViscosity/20;
             elseif strcmp(obj.geometry , 'Square')== 1
                 oilConductance = (oilCrossSectionArea / obj.area)*0.5623 * obj.area^2 *...
                     obj.shapeFactor /network.oilViscosity;
             end              
         end
||||||| merged common ancestors
        function cal_water_conductance(obj, Pc, sig_ow, water_viscosity)
            curvatureRadius = sig_ow / Pc;
            b = zeros(1,4); % Meniscus-Apex distance
            if strcmp(obj.geometry , 'Circle')== 1
                obj.water_area = 0;
                obj.water_conductance = 0;
            else
                dimenlessAreaCorner = zeros(1,4);
                cornerShapeFactor = zeros(1,4);
                scaledDimenlessConducCorner = zeros(1,4);
                dimenlessConducCorner = zeros(1,4);
                halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                for ii = 1:4
                    if ~isnan(halfAngles(ii))
                        b(ii) = (curvatureRadius * cos(obj.receedingContactAngle + halfAngles(ii)) /...
                            sin(halfAngles(ii)));
                        
                        if (halfAngles(ii) + obj.receedingContactAngle) == pi/2
                            dimenlessAreaCorner(ii) = sin(halfAngles(ii))*cos(halfAngles(ii));
                            cornerShapeFactor(ii) = dimenlessAreaCorner(ii) /...
                                (4 * (1 + sin(halfAngles(ii)))^2);
                        else
                            dimenlessAreaCorner(ii) = (sin(halfAngles(ii)) /...
                                cos(obj.receedingContactAngle + halfAngles(ii)))^2 *...
                                ((cos(obj.receedingContactAngle)*...
                                cos(obj.receedingContactAngle + halfAngles(ii)) /...
                                sin(halfAngles(ii)))+ obj.receedingContactAngle +...
                                halfAngles(ii) - pi/2);
                            cornerShapeFactor(ii) = dimenlessAreaCorner(ii) /...
                                (4 * (1 - (sin(halfAngles(ii)) / cos(halfAngles(ii)+...
                                obj.receedingContactAngle))*...
                                (obj.receedingContactAngle + halfAngles(ii) - pi/2))^2);
                        end
                        
                        scaledDimenlessConducCorner(ii) = -15.1794 * cornerShapeFactor(ii)^2 +...
                            7.6307 * cornerShapeFactor(ii) - 0.53488;
                        
                        dimenlessConducCorner(ii) = nthroot(exp((scaledDimenlessConducCorner(ii) +...
                            0.02 * sin(halfAngles(ii) - pi/6)), (1/4/pi - cornerShapeFactor(ii))^(7/8)...
                            / cos(halfAngles(ii) - pi/6)^0.5) * dimenlessAreaCorner(ii)^2);
                    end
                end
                obj.water_area = sum(b.^2 .* dimenlessAreaCorner);
                obj.water_conductance = sum(2 * b.^4 .* dimenlessConducCorner / water_viscosity);
            end
        end
    
        function cal_oil_conductance (obj, oil_viscosity)            
            obj.oil_area = obj.area - obj.water_area;
            if strcmp(obj.geometry , 'Circle')== 1
               obj.oil_conductance = (obj.oil_area / obj.area)*0.5 * obj.area^2 * obj.shapeFactor /oil_viscosity;
            elseif strcmp(obj.geometry , 'Triangle')== 1
               obj.oil_conductance = (obj.oil_area) *3*obj.radius^2/20/oil_viscosity;
            elseif strcmp(obj.geometry , 'Square')== 1
               obj.oil_conductance = (obj.oil_area / obj.area)*0.5623 * obj.area^2 * obj.shapeFactor /oil_viscosity;
            end
        end 

=======
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
>>>>>>> 61a63595e6d7105ef83e8af42fce616fbb52602f
    end
||||||| ad1e304... Change conductance formula
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
         function [waterCrossSectionArea, waterConductance] = calculateWaterConductance(obj, network) 
             halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
             cornerArea = zeros(1,4);
             cornerConductance = zeros(1,4);
             %Raduis of Curvature
             Rc = network.waterViscosity / obj.thresholdPressure;
             for jj = 1:4
                 if ~isnan(halfAngles(jj))
                     % Based on Piri_2005: eq A4 & A5 
                     if (halfAngles(jj) + obj.receedingContactAngle) == pi/2
                          cornerArea(jj) = (Rc * cos(obj.receedingContactAngle + halfAngles(jj))/ sin (halfAngles(jj)))^2 *...
                              sin(halfAngles(jj)) * cos(halfAngles(jj));
                      else
                          cornerArea(jj) = Rc ^2 * (cos(obj.receedingContactAngle)*...
                              (cot(halfAngles(jj)) * cos(obj.receedingContactAngle) - sin(obj.receedingContactAngle))+ ...
                              obj.receedingContactAngle + halfAngles(jj) - pi/2);
                     end   
                     % Based on Piri_2005: eq B(10 - 15)
                     f = 1;
                     F1 = pi/2 - halfAngles(jj) - obj.receedingContactAngle;
                     F2 = cot(halfAngles(jj)) * cos(obj.receedingContactAngle) - sin(obj.receedingContactAngle); 
                     F3 = (pi/2 - halfAngles(jj)) * tan(halfAngles(jj));
                     
                     if (obj.receedingContactAngle <= pi/2 - halfAngles(jj))                         
                         cornerConductance(jj) = (cornerArea(jj)^2 * (1 - sin(halfAngles(jj)))^2 * ...
                             (F2 * cos(obj.receedingContactAngle) - F1) * F3 ^ 2) / ...
                             (12 * network.waterViscosity * ((sin(halfAngles(jj))) * ...
                             (1 - F3) * (F2 + f * F1))^ 2);
                     elseif (obj.receedingContactAngle > pi/2 - halfAngles(jj))
                         cornerConductance(jj) = (cornerArea(jj)^2 * tan(halfAngles(jj))* ...
                             (1 - sin(halfAngles(jj)))^2 * F3 ^ 2) / ...
                             (12 * network.waterViscosity *(sin(halfAngles(jj)))^2*(1 - F3) * (1 + f * F3)^ 2);
%                      else
%                           cornerConductance(jj) = (A1^3 * (1 - sin(halfAngles(jj)))^2 * ...
%                              tan(halfAngles(jj)) * F3 ^ 2) / ...
%                               (12 * network.waterViscosity * (sin(halfAngles(jj)))^2 * ...
%                               (1 - F3) * (1 + f1*F3 - (1- f2*F3) * sqrt(A2/cornerArea(jj)))^ 2);
                     end
                 end
                 waterCrossSectionArea = sum(cornerArea);
                 waterConductance = sum(cornerConductance);
             end
         end
           
         function [oilCrossSectionArea, oilConductance] = calculateOilConductance(obj, network)    
             oilCrossSectionArea = obj.area - obj.waterCrossSectionArea;
             if strcmp(obj.geometry , 'Circle')== 1
                 oilConductance = (oilCrossSectionArea / obj.area)*0.5 * obj.area^2 * ...
                     obj.shapeFactor /network.oilViscosity;
             elseif strcmp(obj.geometry , 'Triangle')== 1
                 oilConductance = oilCrossSectionArea * 3 * obj.radius^2 / network.oilViscosity/20;
             elseif strcmp(obj.geometry , 'Square')== 1
                 oilConductance = (oilCrossSectionArea / obj.area)*0.5623 * obj.area^2 *...
                     obj.shapeFactor /network.oilViscosity;
             end              
         end
    end
=======
        function cal_water_conductance(obj, Pc, sig_ow, water_viscosity)
            curvatureRadius = sig_ow / Pc;
            b = zeros(1,4); % Meniscus-Apex distance
            if strcmp(obj.geometry , 'Circle')== 1
                obj.water_area = 0;
                obj.water_conductance = 0;
            else
                dimenlessAreaCorner = zeros(1,4);
                cornerShapeFactor = zeros(1,4);
                scaledDimenlessConducCorner = zeros(1,4);
                dimenlessConducCorner = zeros(1,4);
                halfAngles = [obj.halfAngle1, obj.halfAngle2,obj.halfAngle3, obj.halfAngle4];
                for ii = 1:4
                    if ~isnan(halfAngles(ii))
                        b(ii) = (curvatureRadius * cos(obj.receedingContactAngle + halfAngles(ii)) /...
                            sin(halfAngles(ii)));
                        
                        if (halfAngles(ii) + obj.receedingContactAngle) == pi/2
                            dimenlessAreaCorner(ii) = sin(halfAngles(ii))*cos(halfAngles(ii));
                            cornerShapeFactor(ii) = dimenlessAreaCorner(ii) /...
                                (4 * (1 + sin(halfAngles(ii)))^2);
                        else
                            dimenlessAreaCorner(ii) = (sin(halfAngles(ii)) /...
                                cos(obj.receedingContactAngle + halfAngles(ii)))^2 *...
                                ((cos(obj.receedingContactAngle)*...
                                cos(obj.receedingContactAngle + halfAngles(ii)) /...
                                sin(halfAngles(ii)))+ obj.receedingContactAngle +...
                                halfAngles(ii) - pi/2);
                            cornerShapeFactor(ii) = dimenlessAreaCorner(ii) /...
                                (4 * (1 - (sin(halfAngles(ii)) / cos(halfAngles(ii)+...
                                obj.receedingContactAngle))*...
                                (obj.receedingContactAngle + halfAngles(ii) - pi/2))^2);
                        end
                        
                        scaledDimenlessConducCorner(ii) = -15.1794 * cornerShapeFactor(ii)^2 +...
                            7.6307 * cornerShapeFactor(ii) - 0.53488;
                        
                        dimenlessConducCorner(ii) = nthroot(exp((scaledDimenlessConducCorner(ii) +...
                            0.02 * sin(halfAngles(ii) - pi/6)), (1/4/pi - cornerShapeFactor(ii))^(7/8)...
                            / cos(halfAngles(ii) - pi/6)^0.5) * dimenlessAreaCorner(ii)^2);
                    end
                end
                obj.water_area = sum(b.^2 .* dimenlessAreaCorner);
                obj.water_conductance = sum(2 * b.^4 .* dimenlessConducCorner / water_viscosity);
            end
        end
    
        function cal_oil_conductance (obj, oil_viscosity)            
            obj.oil_area = obj.area - obj.water_area;
            if strcmp(obj.geometry , 'Circle')== 1
               obj.oil_conductance = (obj.oil_area / obj.area)*0.5 * obj.area^2 * obj.shapeFactor /oil_viscosity;
            elseif strcmp(obj.geometry , 'Triangle')== 1
               obj.oil_conductance = (obj.oil_area) *3*obj.radius^2/20/oil_viscosity;
            elseif strcmp(obj.geometry , 'Square')== 1
               obj.oil_conductance = (obj.oil_area / obj.area)*0.5623 * obj.area^2 * obj.shapeFactor /oil_viscosity;
            end
        end 
>>>>>>> parent of ad1e304... Change conductance formula

    end
end

