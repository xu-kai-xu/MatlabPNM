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
        oil_area
        water_area
        
        %It's better to define the phase conductances and the area of each
        %fluid as a structure for this each element
        fluidConductances
        fluidCrossSectionAreas
        
        receedingContactAngle=20;
        advandingContactAngle
     
    end
    
    methods
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

    end
end

