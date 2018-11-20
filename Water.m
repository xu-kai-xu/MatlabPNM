classdef Water < Fluids
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        viscosity = 0.001;
        sig_ow = 0.02;
        
        
        % Calculated properties
        waterArea;
        crossSectionArea;        
        conductance;
    end
    
    methods
        function obj = Water ()
            obj = obj@Fluids('Water');
        end

    end
    
end

