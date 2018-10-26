classdef Fluids < handle
    %Fluids: contains information related to oil, gas and water 
    
    %   send the conductivity and cross section area of each fluid to each
    %   element and also viscosuty, thermal conductivity, diffusivity of
    %   each phase
    
    properties
        waterViscosity
        oilViscosity
        gasViscosity
        
        waterCrossSectionArea
        oilCrossSectionArea
        gasCrossSectionArea
        waterConductance
        oilConductance
        gasConductance
    end
    
    methods
        function obj = Fluids()
            %Fluids Construct an instance of class Fluids
            %   This class contains fluids properties
            obj.waterViscosity = 0.001;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

