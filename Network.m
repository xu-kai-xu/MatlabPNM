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
        Pc_drain_curve
        
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
        %% Pressure distribution calculation        
        function pressureDistribution (obj, inletPressure, outletPressure, fluid)
            % pressureDistribution Summary of this method goes here
            %   Detailed explanation goes here
            Factor = zeros(obj.numberOfNodes, obj.numberOfNodes);
            B = zeros(obj.numberOfNodes, 1);  
            linkConductance = 0;
            node1Conductance = 0;
            node2Conductance = 0;
            for ii = 1:obj.numberOfLinks
                
                node1Index = obj.Links{ii}.pore1Index;
                node2Index = obj.Links{ii}.pore2Index;      

                % if the link is connected to inlet (index of node 1 is -1 which does not exist) 
                if obj.Links{ii}.isInlet
                    if strcmp(fluid , 'Single_phase')== 1                    
                        linkConductance = obj.Links{ii}.conductance;                             
                        node2Conductance = obj.Nodes{node2Index}.conductance; 
                    elseif strcmp(fluid , 'Two_phase_water')== 1
                        linkConductance = obj.Links{ii}.waterConductance;                              
                        node2Conductance = obj.Nodes{node2Index}.waterConductance;
                    elseif strcmp(fluid , 'Two_phase_oil')== 1                
                        linkConductance = obj.Links{ii}.oilConductance;                   
                        node2Conductance = obj.Nodes{node2Index}.oilConductance;
                    end
                
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        linkConductance) +...
                        0.5 *...
                        ((obj.Links{ii}.pore2Length / node2Conductance)))^-1;
                    
                    Factor(node2Index, node2Index) = Factor(node2Index, node2Index) + nodeLinkSystemConductance;
                    B(node2Index) = nodeLinkSystemConductance * inletPressure;
                % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                elseif obj.Links{ii}.isOutlet
                    if strcmp(fluid , 'Single_phase')== 1                    
                        linkConductance = obj.Links{ii}.conductance;                        
                        node1Conductance = obj.Nodes{node1Index}.conductance; 
                    elseif strcmp(fluid , 'Two_phase_water')== 1
                        linkConductance = obj.Links{ii}.waterConductance;                    
                        node1Conductance = obj.Nodes{node1Index}.waterConductance;
                    elseif strcmp(fluid , 'Two_phase_oil')== 1                
                        linkConductance = obj.Links{ii}.oilConductance;
                        node1Conductance = obj.Nodes{node1Index}.oilConductance; 
                    end
                     nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        linkConductance) +...
                        0.5 *...
                        ((obj.Links{ii}.pore1Length /node1Conductance)))^-1;
                    Factor(node1Index, node1Index) = Factor(node1Index, node1Index) + nodeLinkSystemConductance;
                    B(node1Index) = nodeLinkSystemConductance * outletPressure;
                    
                %if the link is neither inlet nor outlet    
                else
                    if strcmp(fluid , 'Single_phase')== 1                    
                        linkConductance = obj.Links{ii}.conductance;                        
                        node1Conductance = obj.Nodes{node1Index}.conductance;                     
                        node2Conductance = obj.Nodes{node2Index}.conductance; 
                    elseif strcmp(fluid , 'Two_phase_water')== 1
                        linkConductance = obj.Links{ii}.waterConductance;                    
                        node1Conductance = obj.Nodes{node1Index}.waterConductance;                     
                        node2Conductance = obj.Nodes{node2Index}.waterConductance;
                    elseif strcmp(fluid , 'Two_phase_oil')== 1                
                        linkConductance = obj.Links{ii}.oilConductance;
                        node1Conductance = obj.Nodes{node1Index}.oilConductance;                     
                        node2Conductance = obj.Nodes{node2Index}.oilConductance;
                    end
                    nodeLinkSystemConductance = ((obj.Links{ii}.linkLength /...
                        linkConductance) +...
                        0.5 *...
                        ((obj.Links{ii}.pore1Length / node1Conductance) +...
                        (obj.Links{ii}.pore2Length / node2Conductance)))^-1;   
                
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
                obj.Nodes{ii}.Pressure = nodesPressure(ii);      
            end
            
            %assign pressure values to links, since the surface whci
            %flowrate is calculated through might pass through the links
            for ii = 1:obj.numberOfLinks
                if obj.Links{ii}.isInlet
                    obj.Links{ii}.Pressure =...
                        obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure;
                elseif obj.Links{ii}.isOutlet
                    obj.Links{ii}.Pressure =...
                        obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure;                    
                else
                    obj.Links{ii}.Pressure =...
                        (obj.Nodes{obj.Links{ii}.pore1Index}.waterPressure + ...
                        obj.Nodes{obj.Links{ii}.pore2Index}.waterPressure) / 2;
                end
            end
        end
        
        %% Conductance Calculation
         function calculateConductance(obj, Pc)           
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.occupancy == 1 %Node was occupied by oil
                    [obj.Nodes{i}.waterCrossSectionArea, obj.Nodes{i}.waterConductance] =...
                        obj.Nodes{i}.calculateWaterConductance(obj, Pc); 
                    [obj.Nodes{i}.oilCrossSectionArea, obj.Nodes{i}.oilConductance] = ...
                        obj.Nodes{i}.calculateOilConductance(obj);    
                else
                    obj.Nodes{i}.waterCrossSectionArea = obj.Nodes{i}.area;
                    obj.Nodes{i}.waterConductance = obj.Nodes{i}.conductance; 
                    obj.Nodes{i}.oilCrossSectionArea = 0;
                    obj.Nodes{i}.oilConductance = 0;  
                end
            end
            for i = 1:obj.numberOfLinks
                if obj.Links{i}.occupancy == 1 %Link was occupied by oil
                    [obj.Links{i}.waterCrossSectionArea, obj.Links{i}.waterConductance] =...
                        obj.Links{i}.calculateWaterConductance(obj, Pc); 
                    [obj.Links{i}.oilCrossSectionArea, obj.Links{i}.oilConductance] = ...
                        obj.Links{i}.calculateOilConductance(obj);  
                else
                    obj.Links{i}.waterCrossSectionArea = obj.Links{i}.area;
                    obj.Links{i}.waterConductance = obj.Links{i}.conductance; 
                    obj.Links{i}.oilCrossSectionArea = 0;
                    obj.Links{i}.oilConductance = 0;
                end
            end            
        end
        %% Saturation and Conductance Calculation
        function Sw_drain = calculateSaturations(obj, Pc)
            calculateConductance(obj, Pc);
            % Water Saturation Calculation
            waterVolume = 0;            
            for i = 1:obj.numberOfNodes
                if obj.Nodes{i}.occupancy == 1
                     waterVolume = waterVolume + (obj.Nodes{i}.waterCrossSectionArea / obj.Nodes{i}.area) *...
                         obj.Nodes{i}.volume + obj.Nodes{i}.clayVolume;
                else
                     waterVolume = waterVolume + obj.Nodes{i}.volume + obj.Nodes{i}.clayVolume;
                end
            end
            for i = 1:obj.numberOfLinks
                if obj.Links{i}.occupancy == 1
                     waterVolume = waterVolume+ (obj.Links{i}.waterCrossSectionArea / obj.Links{i}.area) * ...
                         obj.Links{i}.volume + obj.Links{i}.clayVolume;
                else
                     waterVolume = waterVolume + obj.Links{i}.volume + obj.Links{i}.clayVolume;
                end  
            end  
            Sw_drain = waterVolume / obj.poreVolume;            
        end
        %% This function calculates the flow rate for each phase in the netwrok
        function calculateFlowRate(obj)
            % location of the surface whcih the flow rate should be
            % calculated through is the half distance of the network
            surfaceLocation = obj.xDimension *0.5;
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
                        
                        %calculate the conductivity of the linkNode system
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
            obj.pressureDistribution(1,0, 'Single_phase');
            obj.calculateFlowRate();
            % unit conversion from m2 to Darcy
            unitConvertor = 1.01325E+12;
            % for pressure difference in the formula the corresponding
            % pressure drop between the vertical surfaces should be
            % calculated (based on Piri B1 formula)
            obj.absolutePermeability = unitConvertor * obj.totalFlowRate * obj.xDimension / (obj.yDimension* obj.zDimension); %/ ()
        end
        
        %% Primary Drainage 
        function obj = PrimaryDrainage(obj) 
            
            % Network is initialy water saturated
            for i = 1:obj.numberOfNodes
                obj.Nodes{i}.occupancy = 0;
            end
            for i = 1:obj.numberOfLinks
                obj.Links{i}.occupancy = 0;
            end            
            
             %% determining the capillary pressure level interval
             Pc_threshold = zeros(obj.numberOfNodes + obj.numberOfLinks,1);
             for i = 1: obj.numberOfNodes                
                 Pc_threshold(i) = obj.Nodes{i}.thresholdPressure;
             end
             for i = obj.numberOfNodes+1 :obj.numberOfLinks                
                 Pc_threshold(i) = obj.Links{i}.thresholdPressure;
             end
             % Pc_interval
             max_Pc = max(Pc_threshold);
             min_Pc = min(Pc_threshold);
             Pc_interval = (max_Pc - min_Pc)/100;
             Pc_drain_max = max_Pc;
             simTimes = Pc_drain_max / Pc_interval;
             fprintf('\nPc_interval is: %f \n', Pc_interval);
             fprintf('Pc_drain_max is: %f \n', Pc_drain_max);
             fprintf('simTimes is: %f \n', simTimes);                
             
             %%     
             % The first row of fill is node/link index, followed by 1 for Nodes and 0 for Links.          
             fill = zeros(obj.numberOfNodes + obj.numberOfLinks,2); % A list for sequense of Node or Link fillings            
             
             %% Updating element saturations and conductances
             a = 0;
             Pc = 1;
             index = 0;
             t = 1;
             obj.Sw_drain(t,1) = calculateSaturations(obj, Pc);  
             %%
             while Pc < Pc_drain_max
                 %% Filling inletNodes by oil          
                 for i = 1:obj.numberOfNodes
                     if obj.Nodes{i}.isInlet                          
                         obj.Nodes{i}.occupancy = 1;
                         a = a + 1;
                         fill(a,1:2) = [obj.Nodes{i}.index, 1];                   
                     end
                 end
                 %% Link with filled Nodes
                 newLinks = zeros(obj.numberOfLinks,2);
                 threshold = zeros(obj.numberOfLinks,1);
                 new =0;
                 for i = 1: obj.numberOfLinks
                     node1Index = obj.Links{i}.pore1Index;
                     node2Index = obj.Links{i}.pore2Index;                     
                     % if the link is connected to inlet (index of node 1 is 0 which does not exist)
                     if obj.Links{i}.isInlet
                         if node1Index == -1 || node1Index == 0
                             index = node2Index;
                         elseif node2Index == -1 || node2Index == 0
                             index = node1Index;
                         end
                         new = new+1;
                         newLinks(new,1) = obj.Links{i}.index;
                         newLinks(new,2) = 1;     
                         threshold(new) =  obj.Links{i}.thresholdPressure;
                         % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                     elseif obj.Links{i}.isOutlet
                         if node1Index == -1 || node1Index == 0
                             index = node2Index;
                         elseif node2Index == -1 || node2Index == 0
                             index = node1Index;
                         end
                         if obj.Nodes{index}.occupancy == 1
                             new = new+1;                                 
                             newLinks(new,1) = obj.Links{i}.index;
                             newLinks(new,2) = 1;     
                             threshold(new) =  obj.Links{i}.thresholdPressure;
                         end
                     elseif obj.Nodes{node1Index}.occupancy ==1 || obj.Nodes{node2Index}.occupancy ==1
                         new = new+1;
                         newLinks(new,1) = obj.Links{i}.index;                         
                         newLinks(new,2) = 1;                 
                         threshold(new) =  obj.Links{i}.thresholdPressure;
                     end                     
                 end 
                 threshold_Throat = zeros(new+10,1);
                 threshold_Throat(1:new,1) = threshold(1:new,1);
                 s=1;
                 while (new > 0)
                     [a,IX]=sort(threshold_Throat(s:new));
                     newLinks = newLinks(IX);
                     node1Index = obj.Links{newLinks(1,1)}.pore1Index;
                     node2Index = obj.Links{newLinks(1,1)}.pore2Index;                     
                     if obj.Links{newLinks(1,1)}.thresholdPressure <=Pc 
                         % if the link is connected to inlet (index of node 1 is 0 which does not exist)
                         if obj.Links{newLinks(1,1)}.isInlet
                             obj.Links{newLinks(1,1)}.occupancy = 1;
                             new = new-1;
                             if node1Index == -1 || node1Index == 0
                                 index = node2Index;
                             elseif node2Index == -1 || node2Index == 0
                                 index = node1Index;
                             end
                             obj.Nodes{index}.occupancy = 1;                             
                             for i=obj.Nodes{index}.connectionNumber
                                 obj.Nodes{index}.connectedLinks(i);
                                 threshold(new+1)= obj.Links{obj.Nodes{index}.connectedLinks(i)}.thresholdPressure;
                                 new = new+1;
                             end
                             % if the link is connected to outlet (index of node 2 is 0 which does not exist)
                         elseif obj.Links{newLinks(1,1)}.isOutlet
                             if node1Index == -1 || node1Index == 0
                                 index = node2Index;
                             elseif node2Index == -1 || node2Index == 0
                                 index = node1Index;
                             end
                             if obj.Nodes{index}.occupancy == 1
                                 obj.Links{newLinks(1,1)}.occupancy = 1;
                                 new = new-1;
                             end
                         elseif(obj.Nodes{node1Index}.occupancy ==1 || ...
                                 obj.Nodes{node2Index}.occupancy ==1)
                             obj.Links{newLinks(1,1)}.occupancy = 1;
                             new = new-1;
                             if obj.Nodes{node1Index}.occupancy ==0
                                 obj.Nodes{node1Index}.occupancy =1;
                                 
                             for i=obj.Nodes{node1Index}.connectionNumber
                                 obj.Nodes{node1Index}.connectedLinks(i);
                                 threshold(new+1)= obj.Links{obj.Nodes{node1Index}.connectedLinks(i)}.thresholdPressure;
                                 new = new+1;
                             end                                
                             elseif obj.Nodes{node2Index}.occupancy ==0
                                 obj.Nodes{node2Index}.occupancy =1;
                                 for i=obj.Nodes{node2Index}.connectionNumber
                                     obj.Nodes{node2Index}.connectedLinks(i);
                                     threshold(new+1)= obj.Links{obj.Nodes{node2Index}.connectedLinks(i)}.thresholdPressure;
                                     new = new+1;
                                 end
                             end
                             s=s+1;
                         end
                         
                     end
                 end
                 
                 %% Updating element saturations and conductances
                 obj.Sw_drain(t,1) = calculateSaturations(obj, Pc);            
                 
                 %% Preparing Pc , Sw & Kr data                     
                 obj.Pc_drain_curve(t,1) = Pc*0.000145037738;                     
                 % Relative Permeability Calculation
                 % [kr_oil(t,1),kr_water(t,1)] = k_rel(1,0);
                 
                 % Pc Step Calculation                     
                 t = t + 1;        
                 %                       if t < length(D)
                 %                           Pc = D(t);
                 %                       else
                 Pc = Pc + Pc_interval;
                 %                       end
                 %                       Pc = Pc + Pc_interval;
             end
             
             plot(obj.Sw_drain,obj.Pc_drain_curve,'--r')
             title('Drainage Cappilary Pressure Curves')
             xlabel('Sw')
             xlim([0 1.05])
             ylabel('Pc (Pa)')
             legend('Drainage Pc')      
        end
        %% CalculateRelativePermeability
        function calculateRelativePermeability(obj)
            %RelativePermeability calculates the relative permeability of
            %oil & water in the network
            %   Detailed explanation goes here
            obj.pressureDistribution(1,0,'Single_phase');
            obj.waterPressure = obj.pressureDistribution(1,0, 'Two_phase_water');
            obj.oilPressure = obj.pressureDistribution(1,0, 'Two_phase_oil');
            obj.calculateFlowRate();
%             % unit conversion from m2 to Darcy
             unitConvertor = 1.01325E+12;
%             % for pressure difference in the formula the corresponding
%             % pressure drop between the vertical surfaces should be
%             % calculated (based on Piri B1 formula)
%             obj.absolutePermeability = unitConvertor * obj.totalFlowRate * obj.xDimension / (obj.yDimension* obj.zDimension); %/ ()
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

