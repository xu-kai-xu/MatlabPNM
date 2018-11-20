classdef Fluids < handle
    %Fluids: contains information related to oil, gas and water 
    
    %   send the conductivity and cross section area of each fluid to each
    %   element and also viscosuty, thermal conductivity, diffusivity of
    %   each phase
    
    properties
        fluidType;
    end
    
    methods (Abstract=true)
        %Conductance = calculateConductance(obj);
        %CrossSectionArea = calculateCrossSectionArea(obj);
    
    end
    
    methods
        function obj = Fluids (Type)
            obj.fluidType = Type;
        end
    end
end

