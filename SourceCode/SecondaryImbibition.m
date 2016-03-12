% clc;clear all
load network.mat
%%
global Pc_interval
%% Calculating throat snap off and pistone like displacement & layer collapse
for i = 1:n_t
    if throat_data(i,5) == 1 % if the throat is oil filled
    [throat_data(i,36) , throat_data(i,37) , throat_data(i,35)] = ...
        throatImbThreshPress(throat_data(i,7),throat_data(i,8:11),throat_data(i,18),...
        throat_data(i,19),throat_data(i,13),throat_data(i,6),sig_ow);
    elseif throat_data(i,5) == 0
        throat_data(i,35) = nan;
    end
    throat_data(i,38:41) = nan; % layer collapse of corners
end
%% Calculating pore snap off & pore body filling pressures & layer collapse
for ii = 1:n_p
    if pore_data(ii,5) == 1 %if the pore is oil filled
    [pore_data(ii,30) , pore_data(ii,31) , pore_data(ii,29)] =...
        poreImbThreshPress(ii,pore_data(ii,7),pore_data(ii,8:11),...
        pore_data(ii,16),pore_data(ii,17),pore_data(ii,12),...
        pore_data(ii,6),sig_ow,pore_data(ii,3));
    elseif pore_data(ii,5) == 0
        pore_data(ii,29) = nan;
    end
    pore_data(ii,32:35) = nan; % layer collapse of corners
end
%% Initial Clustering
[NumberOfClusters, NodeL, LinkL , cluster_A_nums , cluster_B_nums] =...
    clustering(pore_data,throat_data);
%%
A = pore_data(:,30);B = throat_data(:,36);
C = pore_data(:,31);D = throat_data(:,37);
E = zeros(length(A) + length(B) + length(C) + length(D),1);
E(1:length(A)) = A;E(length(A)+1 : length(A)+length(B)) = B;
E(length(A) + length(B) + 1:length(A) + length(B) + length(C)) = C;...
E(length(A) + length(B) + length(C) + 1:end) = D; 
[Nelements,Xcenters] = hist(E,200);
F = Xcenters(Nelements > 5);
F = F(end:-1:1); % Inversing F


Pc_imb = Pc_drain_max;
numOfthroatSnap = 0;
numOfPoreSnap = 0;
numOfPoreFormedLayers = 0;
numOfThroatFormedLayers = 0;
numOfPoreLayerColl = 0;
numOfThroatLayerColl = 0;
numOfThroatPistonLike = 0;
numOfPoreBodyFilling = 0;
t = 1;

while ~isempty(cluster_A_nums) || ~isempty(cluster_B_nums)
    %% invasion check for first layer pores
    for i = 1:n_p
        if pore_data(i,3) == -1 && pore_data(i,31) >= Pc_imb && ...
                (pore_data(i,5) == 1 && ~any( pore_data(i,32:35))) && ...
                any(NodeL(i) == cluster_A_nums(:))
            if ~isnan(pore_data(i,29))
                pore_data(i,5) = 1; % an oil sandwich layer exist
                numOfPoreBodyFilling = numOfPoreBodyFilling + 1;
                % Calculating Pc layer collapse in the pore
                [pore_data(i,32:35),pore_data(i,29)] =...
                    PcLayerCollapse(pore_data(i,7) , ...
                    pore_data(i,8:11),pore_data(i,17),...
                    Pc_imb , pore_data(i,16));
             else
                pore_data(i,5) = 0; % inaved the pore with water
                numOfPoreBodyFilling = numOfPoreBodyFilling + 1;
            end
        end
    end
    %% Percolation List Construction
    [~, NodeL, LinkL , cluster_A_nums , cluster_B_nums] =...
        clustering(pore_data,throat_data);
    
%     [pore_data, throat_data] = clusterTrapProp(Pc_imb ,pore_data, throat_data, NodeL, ...
%         LinkL, cluster_A_nums , cluster_B_nums);

    percList = [];
    a = 1;
    for ii = 1:n_t
        if (any(LinkL(ii) == cluster_A_nums(:)) && ~any(throat_data(ii,38:41))) || ...
                (any(LinkL(ii) == cluster_B_nums(:)) && ~any(throat_data(ii,38:41)))
            if (pore_data(throat_data(ii,21),5) == 0 || (pore_data(throat_data(ii,21),5) == 1 ...
                    && any(pore_data(throat_data(ii,21),32:35)))) ||...
                    (pore_data(throat_data(ii,22),5) == 0 || (pore_data(throat_data(ii,22),5) == 1 ...
                    && any(pore_data(throat_data(ii,22),32:35))))
                
                percList(a,1) = ii;
                percList(a,2) = 2;
                percList(a,3) = 2;
                percList(a,4) = throat_data(ii,37);
                a = a + 1;
            elseif ~isnan(throat_data(ii,36))% if the throat is non circular
                percList(a,1) = ii;
                percList(a,2) = 2;
                percList(a,3) = 1;
                percList(a,4) = throat_data(ii,36);
                a = a + 1;                
            end
        end
    end
    for ii = 1:n_p
        if ((any(NodeL(ii) == cluster_A_nums(:)) && ~any(pore_data(ii,32:35))) ||...
                (any(NodeL(ii) == cluster_B_nums(:)) && ~any(pore_data(ii,32:35))))...
               && pore_data(ii,3) ~= -1
            
            filledThroats = 0;
            for j = 1:length(nonzeros(pore_data(ii,41:end)))
                if throat_data(pore_data(ii,j+40),5) == 0 || ...
                        (any(LinkL(pore_data(ii,j+40)) == cluster_A_nums(:)) &&...
                        any(LinkL(pore_data(ii,j+40)) == cluster_B_nums(:)) &&...
                        any(throat_data(pore_data(ii,j+40),38:41))) 
%                         (throat_data(pore_data(ii,j+40),5) == 1 &&...
%                         any(throat_data(pore_data(ii,j+40),38:41)))
                    filledThroats = filledThroats + 1;
                end
            end
            if filledThroats ~= 0 % yani agar in pore, throat havi ab dasht
                percList(a,1) = ii;
                percList(a,2) = 1;
                percList(a,3) = 2;
                percList(a,4) = pore_data(ii,31); % pore body filling threshold pressure
                a = a + 1;
            elseif ~isnan(pore_data(ii,30)) % if the pore is not circular which means snap off threshold pressure is not nan
                percList(a,1) = ii;
                percList(a,2) = 1;
                percList(a,3) = 1;
                percList(a,4) = pore_data(ii,30); % snap off threshold pressure
                a = a + 1;
            end
        end
    end
    
    %% Percolation Section
    while ~isempty(percList)
        [~, NodeL, LinkL , cluster_A_nums , cluster_B_nums] =...
            clustering(pore_data,throat_data);
        % Descending sorting of threshold pressures
        percList = sortrows(percList,-4);
    
        if percList(1,2) == 2 % if the first element is a throat
            if any(LinkL(percList(1,1)) == cluster_A_nums(:)) || ...
                    any(LinkL(percList(1,1)) == cluster_B_nums(:))
                
                if percList(1,3) == 1 && Pc_imb < percList(1,4)
                    throat_data(percList(1,1),5) = 0;
                    throat_data(percList(1,1),35) = nan;
                    
                    numOfthroatSnap = numOfthroatSnap + 1;

                    % Eliminating snap off feature of the  pores attached to
                    % the water invaded throat from percList
                    eliminatePoreIndex1 = find(percList(:,1) == throat_data(percList(1,1),21) &...
                        percList(:,2) == 1 & percList(:,3) == 1);
                    percList(eliminatePoreIndex1,:) = [];
                    eliminatePoreIndex2 = find(percList(:,1) == throat_data(percList(1,1),22) &...
                        percList(:,2) == 1 & percList(:,3) == 1);
                    percList(eliminatePoreIndex2,:) = [];
                    if ~any(pore_data(throat_data(percList(1,1),21),32:35)) && ... % no layer exist in the pore
                            (any(NodeL(throat_data(percList(1,1),21)) == cluster_A_nums(:)) || ...
                            any(NodeL(throat_data(percList(1,1),21)) == cluster_B_nums(:)))

                        % shayad niaz bashad dar inja shart inke pore marboot
                        % be laye aval nabashad ra ham biavarim

                        % if the pore is not in the percList
                        if isempty(find(percList(:,1) == throat_data(percList(1,1),21) &...
                        percList(:,2) == 1 & percList(:,3) == 2, 1))

                            percList(end+1,:) = [pore_data(throat_data(percList(1,1),21),1) , 1 , 2 ,...
                                pore_data(throat_data(percList(1,1),21),31)];
                        end
                    end
                    % elseif nagozashtam choon shyad har do halat barqarar
                    % bashad !
                    if ~any(pore_data(throat_data(percList(1,1),22),32:35)) && ...
                            (any(NodeL(throat_data(percList(1,1),22)) == cluster_A_nums(:)) || ...
                            any(NodeL(throat_data(percList(1,1),22)) == cluster_B_nums(:)))
                        if isempty(find(percList(:,1) == throat_data(percList(1,1),22) &...
                        percList(:,2) == 1 & percList(:,3) == 2, 1))

                            percList(end+1,:) = [pore_data(throat_data(percList(1,1),22),1) , 1 , 2 ,...
                                pore_data(throat_data(percList(1,1),22),31)];
                        end
                    end
                    percList(1,:) = [];
                elseif percList(1,3) == 2 && Pc_imb < percList(1,4)
                    if ~isnan(throat_data(percList(1,1),35)) % If layer formation is possible
                        throat_data(percList(1,1),5) = 1; % make the throat oil type
                        numOfThroatPistonLike = numOfThroatPistonLike + 1;
                        % calculating Pc collapse of the layers
                        [throat_data(percList(1,1),38:41),throat_data(percList(1,1),35)] =...
                            PcLayerCollapse(throat_data(percList(1,1),7), ...
                            throat_data(percList(1,1),8:11),throat_data(percList(1,1),19),...
                            Pc_imb , throat_data(percList(1,1),18));
                    else
                        throat_data(percList(1,1),5) = 0;
                        numOfThroatPistonLike = numOfThroatPistonLike + 1;
                    end

                    % Eliminating snap off feature of pores attached to the
                    % water invaded throat in percList
                    eliminatePoreIndex1 = find(percList(:,1) == throat_data(percList(1,1),21) &...
                        percList(:,2) == 1 & percList(:,3) == 1);
                    percList(eliminatePoreIndex1,:) = [];
                    eliminatePoreIndex2 = find(percList(:,1) == throat_data(percList(1,1),22) &...
                        percList(:,2) == 1 & percList(:,3) == 1);
                    percList(eliminatePoreIndex2,:) = [];
                    if ~any(pore_data(throat_data(percList(1,1),21),32:35)) && ...
                            (any(NodeL(throat_data(percList(1,1),21)) == cluster_A_nums(:)) || ...
                            any(NodeL(throat_data(percList(1,1),21)) == cluster_B_nums(:)))

                        % Updating pore body filling of the pore
                        [pore_data(throat_data(percList(1,1),21),31) ,...
                            pore_data(throat_data(percList(1,1),21),29)] = poreBodyFillingThreshPress...
                            (pore_data(throat_data(percList(1,1),21),1),pore_data(throat_data(percList(1,1),21),3),...
                            pore_data(throat_data(percList(1,1),21),17),pore_data(throat_data(percList(1,1),21),7),...
                            pore_data(throat_data(percList(1,1),21),6),pore_data(throat_data(percList(1,1),21),12),...
                            pore_data(throat_data(percList(1,1),21),8:11),pore_data(throat_data(percList(1,1),21),16));

                        percList(end+1,:) = [pore_data(throat_data(percList(1,1),21),1) , 1 , 2 ,...
                            pore_data(throat_data(percList(1,1),21),31)];
                    end
                    % elseif nagozashtam choon shyad har do halat barqarar
                    % bashad !
                    if ~any(pore_data(throat_data(percList(1,1),22),32:35)) && ...
                            (any(NodeL(throat_data(percList(1,1),22)) == cluster_A_nums(:)) || ...
                            any(NodeL(throat_data(percList(1,1),22)) == cluster_B_nums(:)))
                        
                        % Updating pore body filling of the pore
                        [pore_data(throat_data(percList(1,1),22),31) , ...
                            pore_data(throat_data(percList(1,1),22),29)] = poreBodyFillingThreshPress...
                            (pore_data(throat_data(percList(1,1),22),1),pore_data(throat_data(percList(1,1),22),3),...
                            pore_data(throat_data(percList(1,1),22),17),pore_data(throat_data(percList(1,1),22),7),...
                            pore_data(throat_data(percList(1,1),22),6),pore_data(throat_data(percList(1,1),22),12),...
                            pore_data(throat_data(percList(1,1),22),8:11),pore_data(throat_data(percList(1,1),22),16));

                        percList(end+1,:) = [pore_data(throat_data(percList(1,1),22),1) , 1 , 2 ,...
                            pore_data(throat_data(percList(1,1),22),31)];
                    end
                    percList(1,:) = [];
                else
                   percList(1,:) = []; 
                end
            else
                percList(1,:) = [];
            end
            
        elseif percList(1,2) == 1 % if the first element is a pore
            if any(NodeL(percList(1,1)) == cluster_A_nums(:)) || ...
                    any(NodeL(percList(1,1)) == cluster_B_nums(:))

                if percList(1,3) == 1 && Pc_imb < percList(1,4)
                    pore_data(percList(1,1),5) = 0;
                    pore_data(percList(1,1),29) = nan;
                    
                    numOfPoreSnap = numOfPoreSnap + 1;

                    for ii = 1:length(nonzeros(pore_data(percList(1,1),41:end)))
                        eliminatedthroat = find(percList(:,1) == throat_data(pore_data(percList(1,1),40+ii),1) &...
                            percList(:,2) == 2 & percList(:,3) == 1);
                        percList(eliminatedthroat,:) = [];

                        if ~any(throat_data(pore_data(percList(1,1),40+ii),38:41)) && ...
                                (any(LinkL(throat_data(pore_data(percList(1,1),40+ii),1)) == cluster_A_nums(:)) || ...
                                any(LinkL(throat_data(pore_data(percList(1,1),40+ii),1)) == cluster_B_nums(:)))

                            percList(end + 1,:) = [throat_data(pore_data(percList(1,1),40+ii),1) , 2 ,2 ,...
                                throat_data(pore_data(percList(1,1),40+ii),37)];
                        end
                    end
                    percList(1,:) = [];

                elseif percList(1,3) == 2 && Pc_imb < percList(1,4)

                    if ~isnan(pore_data(percList(1,1),29)) % If layer formation is impossible
                        pore_data(percList(1,1),5) = 1; % make the pore oil type
                        numOfPoreBodyFilling = numOfPoreBodyFilling + 1;
                        
                        [pore_data(percList(1,1),32:35),pore_data(percList(1,1),29)] =...
                            PcLayerCollapse(pore_data(percList(1,1),7) , ...
                            pore_data(percList(1,1),8:11),pore_data(percList(1,1),17),...
                            Pc_imb , pore_data(percList(1,1),16));
                    else % if no layer will not form
                        pore_data(percList(1,1),5) = 0; % make the pore water type
                        numOfPoreBodyFilling = numOfPoreBodyFilling + 1;
                    end

                    for ii = 1:length(nonzeros(pore_data(percList(1,1),41:end)))
                        eliminatedthroat = find(percList(:,1) == ...
                            throat_data(pore_data(percList(1,1),40+ii),1) & percList(:,2) == 2 &...
                            percList(:,3) == 1);
                        percList(eliminatedthroat,:) = [];

                        if ~any(throat_data(pore_data(percList(1,1),40+ii),38:41)) && ...
                                (any(LinkL(throat_data(pore_data(percList(1,1),40+ii),1)) == cluster_A_nums(:)) || ...
                                any(LinkL(throat_data(pore_data(percList(1,1),40+ii),1)) == cluster_B_nums(:)))
                            percList(end + 1,:) = [throat_data(pore_data(percList(1,1),40+ii),1) , 2 ,2 ,...
                                throat_data(pore_data(percList(1,1),40+ii),37)];
                        end
                    end
                    percList(1,:) = [];
                else
                    percList(1,:) = [];
                end
            else
                percList(1,:) = [];
            end
        end
    end
    %% Updating Pc collapse of the layers
    for ii = 1:n_p
        if ~isnan(pore_data(ii,29)) && any(pore_data(ii,32:35))... %&& pore_data(ii,5) == 1
                && (any(NodeL(ii) == cluster_A_nums(:)) || any(NodeL(ii) == cluster_B_nums(:)))
            % Updating Pc of layer collapse
%             pore_data(ii,32:35) = PcLayerCollapseUpdate(pore_data(ii,32:35),...
%                 pore_data(ii,8:11),pore_data(ii,17),Pc_imb,pore_data(ii,29),pore_data(ii,7));
            %Cheking layer collapse
            for jj = 1:4
                if ~isnan(pore_data(ii,jj + 31)) && pore_data(ii,jj + 31) > Pc_imb
                    pore_data(ii,jj + 31) = nan;
                    numOfPoreLayerColl = numOfPoreLayerColl + 1;

                end
            end
            if ~any(pore_data(ii,32:35))
                pore_data(ii,5) = 0;
                pore_data(ii,29) = nan;
            end
        end
    end
    for ii = 1:n_t
        if ~isnan(throat_data(ii,35)) && any(throat_data(ii,38:41))...% && throat_data(ii,5) == 1
                && (any(LinkL(ii) == cluster_A_nums(:)) || any(LinkL(ii) == cluster_B_nums(:)))
            % Updating Pc of layer collapse
%             throat_data(ii,38:41) = PcLayerCollapseUpdate(throat_data(ii,38:41),...
%                 throat_data(ii,8:11),throat_data(ii,19),Pc_imb,throat_data(ii,35),throat_data(ii,7));
            %Cheking layer collapse
            for jj = 1:4
                if ~isnan(throat_data(ii,37 + jj)) && throat_data(ii,37 + jj) > Pc_imb
                    throat_data(ii,37 + jj) = nan;
                    numOfThroatLayerColl = numOfThroatLayerColl + 1;
                end
            end
            if ~any(throat_data(ii,38:41))
                throat_data(ii,5) = 0;
                throat_data(ii,35) = nan;
            end
        end
    end
    r=1;
    
%     [NumberOfClusters, NodeL, LinkL , cluster_A_nums , cluster_B_nums] =...
%         clustering(pore_data,throat_data);
    %% Updating saturations and conductances
    for ii = 1:n_p
        if pore_data(ii,5) == 0 || (pore_data(ii,5) == 1 && any(pore_data(ii,32:35)))
            % Piri 
            [pore_data(ii,27),pore_data(ii,26),pore_data(ii,22),pore_data(ii,21)] =...
                conduct_imb(Pc_imb,pore_data(ii,16),pore_data(ii,17),pore_data(ii,8:11),...
                pore_data(ii,7),pore_data(ii,14),pore_data(ii,29),pore_data(ii,25),...
                pore_data(ii,6),pore_data(ii,12),pore_data(ii,29),pore_data(ii,32:35));            
            % Patzek
%             [pore_data(ii,27),pore_data(ii,26),pore_data(ii,22),pore_data(ii,21)] =...
%                 conduct_imb_7(Pc_imb,pore_data(ii,16),pore_data(ii,17),pore_data(ii,8:11),...
%                 pore_data(ii,7),pore_data(ii,14),pore_data(ii,29),pore_data(ii,25),...
%                 pore_data(ii,6),pore_data(ii,12),pore_data(ii,29),pore_data(ii,32:35));
%         elseif pore_data(ii,5) == 1 && ~any(pore_data(ii,32:35)) %&& ...
%    %             (any(NodeL(ii) == cluster_A_nums(:)) || any(NodeL(ii) == cluster_B_nums(:)))
%             % Piri
% %             [pore_data(ii,27),pore_data(ii,26),pore_data(ii,22),pore_data(ii,21)] =...
% %                 conduct_imb_6_1(Pc_imb,pore_data(ii,16),pore_data(ii,8:11),...
% %                 pore_data(ii,7),pore_data(ii,6),pore_data(ii,14),sig_ow,...
% %                 pore_data(ii,12), pore_data(ii,17));
%             % Patzek
%             [pore_data(ii,27),pore_data(ii,26),pore_data(ii,22),pore_data(ii,21)] =...
%                 conduct_imb_7_1(Pc_imb,pore_data(ii,16),pore_data(ii,8:11),...
%                 pore_data(ii,7),pore_data(ii,6),pore_data(ii,14),sig_ow,...
%                 pore_data(ii,12));
        end
            % Valvatne
            
            % Update Fluid Volume
            [pore_data(ii,24),pore_data(ii,23)] = Cal_fluid_vol...
                (pore_data(ii,22),pore_data(ii,21),pore_data(ii,13));
    end
    for ii = 1:n_t
        if throat_data(ii,5) == 0 || (throat_data(ii,5) == 1 && any(throat_data(ii,38:41))) % water type
            % Piri
            [throat_data(ii,33),throat_data(ii,32),throat_data(ii,28),throat_data(ii,27)] =...
                conduct_imb(Pc_imb,throat_data(ii,18),throat_data(ii,19),throat_data(ii,8:11),...
                throat_data(ii,7),throat_data(ii,15),throat_data(ii,35),throat_data(ii,31),...
                throat_data(ii,6),throat_data(ii,13),throat_data(ii,35),throat_data(ii,38:41));
            % Patzek
%             [throat_data(ii,33),throat_data(ii,32),throat_data(ii,28),throat_data(ii,27)] =...
%                 conduct_imb_7(Pc_imb,throat_data(ii,18),throat_data(ii,19),throat_data(ii,8:11),...
%                 throat_data(ii,7),throat_data(ii,15),throat_data(ii,35),throat_data(ii,31),...
%                 throat_data(ii,6),throat_data(ii,13),throat_data(ii,35),throat_data(ii,38:41));
%         elseif throat_data(ii,5) == 1 && ~any(throat_data(ii,38:41))% && ...
%     %            (any(LinkL(ii) == cluster_A_nums(:)) || any(LinkL(ii) == cluster_B_nums(:)))% Oil type
%             % Piri
% %             [throat_data(ii,33),throat_data(ii,32),throat_data(ii,28), throat_data(ii,27)] =...
% %                 conduct_imb_6_1(Pc_imb,throat_data(ii,18),throat_data(ii,8:11),throat_data(ii,7),...
% %                 throat_data(ii,6),throat_data(ii,15),sig_ow,throat_data(ii,13), throat_data(ii,19));
%             % Patzek
%             [throat_data(ii,33),throat_data(ii,32),throat_data(ii,28), throat_data(ii,27)] =...
%                 conduct_imb_7_1(Pc_imb,throat_data(ii,18),throat_data(ii,8:11),throat_data(ii,7),...
%                 throat_data(ii,6),throat_data(ii,15),sig_ow,throat_data(ii,13));
        end
            % Valvatne
            
            % Updating Fluid volume
            [throat_data(ii,30),throat_data(ii,29)] = Cal_fluid_vol(throat_data(ii,28),...
                throat_data(ii,27),throat_data(ii,14));
    end
    
    
    %% Relative permeability calculation
    [kr_oil_imb(t,1),kr_water_imb(t,1)] = k_rel_imb(1,0,NodeL,cluster_A_nums);
    
    %% Water saturation calculation
    water_vol = sum(pore_data(:,23)) + sum(throat_data(:,29));
    Sw_imb(t,1) = water_vol / totalPV

    %%
    Pc_imb_curve(t,1) = Pc_imb;
    t = t + 1
%     if t < length(F)
%         Pc_imb = -F(t);
%     else
%         Pc_imb = Pc_imb - Pc_interval;
%     end

%     if t == 2
%         Pc_imb = Pc_drain_max/3;
%     else
%         Pc_imb = Pc_imb - 0.3*Pc_interval;
%     end

   if ~isempty(cluster_A_nums) || ~isempty(cluster_B_nums)
        Pc_imb = F(t - 1);
    else
        Pc_imb = Pc_imb - 1*Pc_interval;
    end
 
    [NumberOfClusters, NodeL, LinkL , cluster_A_nums , cluster_B_nums] =...
        clustering(pore_data,throat_data);    
    
    [pore_data, throat_data] = clusterTrapProp(Pc_imb ,pore_data, throat_data, NodeL, ...
        LinkL, cluster_A_nums , cluster_B_nums);

end
Pc_imb_curve = Pc_imb_curve * 0.000145037738;

%% Ploting Pc & Kr
figure('name','Primary Dranage & Secondary Imbibition Cappilary Pressure & Relative Permeability Curves',...
    'units','normalized','outerposition',[0 0 1 1])
subplot(2,2,[1 3]);
grid on
plot(Sw_drain,Pc_drain_curve,'--r')
hold on
plot(Sw_imb,Pc_imb_curve,'-b')
title('Drainage & Imbibition Cappilary Pressure Curves')
xlabel('Sw')
xlim([0 1.05])
ylabel('Pc (Pa)')
legend('Drainage','Imbibition')

subplot(2,2,2);
plot(Sw_drain,kr_water,'-r',Sw_drain,kr_oil,'--b')
title('Drainage Relative Permeability Curves')
xlabel('Sw')
xlim([0 1.05])
ylabel('Reative Permeability')
ylim([0 1.05])
legend('Water Relative Permeability','Oil Relative Permeability','Location','West')

subplot(2,2,4);
plot(Sw_imb,kr_water_imb,'-r',Sw_imb,kr_oil_imb,'--b')
title('Imbibition Relative Permeability Curves')
xlabel('Sw')
xlim([0 1.05])
ylabel('Reative Permeability')
ylim([0 1.05])
legend('Water Relative Permeability','Oil Relative Permeability')