classdef Network < handle & Fluids
    %Network Summary of this class goes here
    %   This class contain nodes and links of the network
    
    properties
        Nodes
        Links
        xDimension
        yDimension
        zDimension
        numberOfLinks
        numberOfNodes
        
        Porosity
        totalFlowRate
        absolutePermeability
        poreVolume
        Sw_drain
    end
    
    methods
        %% Cunstructor function
        function obj = Network(fileName)
            %Network Construct an instance of this class
            %   Detailed explanation goes here
            
%             fileName = 'Berea';
            
            % Opening the files
            link_1_fileID = fopen(strcat(fileName, '_link1.dat'));
            obj.numberOfLinks = str2num(fgetl(link_1_fileID));
            link_2_fileID = fopen(strcat(fileName, '_link2.dat'));
            
            node_2_fileID = fopen(strcat(fileName, '_node2.dat'));            
            node_1_fileID = fopen(strcat(fileName, '_node1.dat'));
            temp = str2num(fgetl(node_1_fileID));
            obj.numberOfNodes = temp(1);
            
            % Network dimension
            obj.xDimension = temp(2);
            obj.yDimension = temp(3);
            obj.zDimension = temp(4);
            
            % Initializing Nodes and Links parameters
            obj.Nodes = cell(obj.numberOfNodes,1);
            obj.Links = cell(obj.numberOfLinks,1);
            
            % 
            for i = 1:obj.numberOfNodes
                node_1_values = str2num(fgetl(node_1_fileID));
                node_2_values = str2num(fgetl(node_2_fileID));
                obj.Nodes{i} = Node(node_1_values(1),... %pore index
                                    node_1_values(2),... % pore x coordinate
                                    node_1_values(3),... % pore y coordinate
                                    node_1_values(4),... % pore z coordinate
                                    node_1_values(5),... %pore connection number
                                    node_1_values(6:end),... % inlet-outlet status and connected link index
                                    node_2_values(2),... % pore volume
                                    node_2_values(3),... % pore radius  
                                    node_2_values(4),... % pore shape factor 
                                    node_2_values(5)); % pore clay volume               
            end        
            
            for i = 1:obj.numberOfLinks
               link_1_values = str2num(fgetl(link_1_fileID));
               link_2_values = str2num(fgetl(link_2_fileID));
               obj.Links{i} = Link(link_1_values(1),... %index 
                                    link_1_values(2),... %pore1Index,... 
                                    link_1_values(3),... %pore2Index,...
                                    link_1_values(4),... %radius,...
                                    link_1_values(5),... %shapeFactor,...
                                    link_1_values(6),... %length,...
                                    link_2_values(4),... %pore1Length,...
                                    link_2_values(5),... %pore2Length,...
                                    link_2_values(6),... %linkLength,...
                                    link_2_values(7),... %volume,...
                                    link_2_values(8)); %clayVolume
            end
            
            %closing the files
            fclose(link_1_fileID); fclose(link_2_fileID);
            fclose(node_1_fileID); fclose(node_2_fileID);    
            
        end
        %% Porosity calculation
        function obj = calculatePorosity(obj)
            nodesVolume = 0;
            linksVolume = 0;
            for ii = 1:obj.numberOfNodes
                nodesVolume = nodesVolume + (obj.Nodes{ii}.volume);
            end
            for ii = 1:obj.numberOfLinks
                linksVolume = linksVolume + (obj.Links{ii}.volume); 
            end
            obj.Porosity = (linksVolume + nodesVolume) / (obj.xDimension * obj.yDimension * obj.zDimension);  
            obj.poreVolume = linksVolume + nodesVolume;
        end
        %% Conductance Calculation
         function calculateConductance(obj)           
            for ii = 1:obj.numberOfNodes                
             [obj.Nodes{ii}.waterArea, obj.Nodes{ii}.waterConductance] =...
                 obj.Nodes{ii}.calculateWaterConductance(obj); 
             [obj.Nodes{ii}.oilArea, obj.Nodes{ii}.oilConductance] = ...
                 obj.Nodes{ii}.calculateOilConductance(obj); 
            end
            for ii = 1:obj.numberOfLinks                
             [obj.Links{ii}.waterArea, obj.Links{ii}.waterConductance] = ...
                 obj.Links{ii}.calculateWaterConductance(obj); 
%              [obj.Links{ii}.oilArea, obj.Links{ii}.oilConductance] = ...
%                  obj.Links{ii}.calculateOilConductance(obj);
            end                    
        end
        %% Saturation Calculation
        function calculateSaturations(obj)
            % Water Saturation Calculation
            waterVolume = 0;                     
            for ii = 1:obj.numberOfLinks                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;
                 % if the link is connected to inlet (index of node 1 is -1 which does not exist) 
                if obj.Links{ii}.isInlet
                waterVolume = waterVolume+ obj.Links{ii}.waterArea * obj.Links{ii}.linkLength + ...                   
                    obj.Nodes{node2Index}.waterArea* obj.Links{ii}.pore2Length;
%                 totalPV = totalPV+ obj.Links{ii}.volume + obj.Nodes{node2Index}.volume;
                 % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                elseif obj.Links{ii}.isOutlet
                 %if the link is neither inlet nor outlet  
                 waterVolume = waterVolume+ obj.Links{ii}.waterArea * obj.Links{ii}.linkLength + ...
                     obj.Nodes{node1Index}.waterArea* obj.Links{ii}.pore1Length;
%                  totalPV = totalPV + obj.Links{ii}.volume + obj.Nodes{node1Index}.volume;
                else
                    waterVolume = waterVolume+ obj.Links{ii}.waterArea * obj.Links{ii}.linkLength + ...
                        obj.Nodes{node1Index}.waterArea* obj.Links{ii}.pore1Length + ...
                        obj.Nodes{node2Index}.waterArea* obj.Links{ii}.pore2Length;
%                     totalPV =  totalPV+ obj.Links{ii}.volume + ...
%                         obj.Nodes{node1Index}.volume + obj.Nodes{node2Index}.volume;
                end
                    
            end
            obj.Sw_drain = waterVolume / obj.poreVolume;
            
        end
        %% Pressure distribution calculation        
        function pressureDistribution (obj, inletPressure, outletPressure)
            % pressureDistribution Summary of this method goes here
            %   Detailed explanation goes here
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            B = zeros(obj.numberOfNodes, 1);  
            
            for ii = 1:obj.numberOfLinks
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;

                % if the link is connected to inlet (index of node 1 is -1 which does not exist) 
                if obj.Links{ii}.isInlet
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance) +...
                        0.5 *...
                        ((obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance)))^-1;
                    
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + nodeLinkSystemConductance;
                    B(node2Index) = nodeLinkSystemConductance * inletPressure;
                % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                elseif obj.Links{ii}.isOutlet
                     nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance) +...
                        0.5 *...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance)))^-1;
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + nodeLinkSystemConductance;
                    B(node1Index) = nodeLinkSystemConductance * outletPressure;
                    
                %if the link is neither inlet nor outlet    
                else
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        obj.Links{ii}.conductance) +...
                        0.5 *...
                        ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance) +...
                        (obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance)))^-1;   
                
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + nodeLinkSystemConductance;
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + nodeLinkSystemConductance;
                    Factor(node1Index, node2Index) = Factor(node1Index, node2Index) - nodeLinkSystemConductance;
                    Factor(node2Index, node1Index) = Factor(node2Index, node1Index) - nodeLinkSystemConductance;
                   
                end     
            end
            
            % using Preconditioned conjugate gradients method to solve the
            % pressure distribution 
            nodesPressure = pcg(Factor,B,1e-3,300);
            
            %assign the pressure values to each node
            for ii = 1:obj.numberOfNodes
                obj.Nodes{ii}.waterPressure = nodesPressure(ii);      
            end
            
            %assign pressure values to links, since the surface whci
            %flowrate is calculated through might pass through the links
            for ii = 1:obj.numberOfLinks
                if obj.Links{ii}.isInlet
                    obj.Links{ii}.waterPressure =...
                        obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure;
                elseif obj.Links{ii}.isOutlet
                    obj.Links{ii}.waterPressure =...
                        obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure;                    
                else
                    obj.Links{ii}.waterPressure =...
                        (obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure + ...
                        obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure) / 2;
                end
            end
        end
        %% This function calculates the flow rate for each phase in the netwrok
        function calculateFlowRate(obj)
            % location of the surface whcih the flow rate should be
            % calculated through is the half distance of the network
            surfaceLocation = obj.xDimension / 2;
            flowRate = 0;
            
            %search through all the links
            for ii = 1:obj.numberOfLinks
                %the link should not be nether inlet nor outlet becasue it
                %makes problem in the index of the connected nodes
                if ~obj.Links{ii}.isInlet && ~obj.Links{ii}.isOutlet
                    
                    node1Index = obj.Links{ii}.pore1Index;
                    node2Index = obj.Links{ii}.pore2Index;
                   %if the two connected nodes pass through the sufrace
                    %count the flow of fluid passing the link connecting
                    %them
                    if obj.Nodes{node1Index}.x_coordinate < surfaceLocation && ...
                            obj.Nodes{node2Index}.x_coordinate > surfaceLocation   
                        
                        %calculate the conductivity of the linNode system
                        nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                            obj.Links{ii}.conductance) +...
                            0.5 *...
                            ((obj.Links{ii}.pore1Length / obj.Nodes{node1Index}.conductance) +...
                            (obj.Links{ii}.pore2Length / obj.Nodes{node2Index}.conductance)))^-1;
                        
                        % calculate the flow rate of the fluid
                        flowRate = flowRate + ...
                            nodeLinkSystemConductance * ...
                            (obj.Nodes{node1Index}.waterPressure - ...
                            obj.Nodes{node2Index}.waterPressure);
                        
                    end
                end 
            end
            obj.totalFlowRate = flowRate;           
        end
        %% Calculate AbsolutePermeability 
        function calculateAbsolutePermeability(obj)
            %AbsolutePermeability calculates the absolute permeability of
            %the network
            %   Detailed explanation goes here
            obj.pressureDistribution(1,0);
            obj.calculateFlowRate();
            % unit conversion from m2 to Darcy
            unitConvertor = 1.01325E+12;
            % for pressure difference in the formula the corresponding
            % pressure drop between the vertical surfaces should be
            % calculated (based on Piri B1 formula)
            obj.absolutePermeability = unitConvertor * 0.001 * obj.totalFlowRate * obj.xDimension / (obj.yDimension* obj.zDimension); %/ ()
        end     
        %% vtk file generation
        function vtkOutput(obj)
            vtkFileID = fopen('output.vtk','w');
            if vtkFileID == -1
                error('Cannot open file for writing.');
            end
            title = 'output';
            fprintf ( vtkFileID, '# vtk DataFile Version 2.0\n' );
            fprintf ( vtkFileID, '%s\n', title );
            fprintf ( vtkFileID, 'ASCII\n' );
            fprintf ( vtkFileID, '\n' );
            fprintf ( vtkFileID, 'DATASET POLYDATA\n' );
            fprintf ( vtkFileID, 'POINTS %d double\n', obj.numberOfNodes );
            for i = 1:obj.numberOfNodes
                fprintf( vtkFileID,'%d %d %d \n', obj.Nodes{i}.x_coordinate, obj.Nodes{i}.y_coordinate, obj.Nodes{i}.z_coordinate );
            end
            
        end
    end
end

