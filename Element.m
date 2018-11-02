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
        oil_area
        water_area
        
        %It's better to define the phase conductances and the area of each
        %fluid as a structure for this each element
        waterConductance
        oilConductance
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
        
        
    
        function obj = cal_oil_conductance (obj, oil_viscosity)            
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

