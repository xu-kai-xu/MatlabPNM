% clc; clear all
load network.mat
global pore_data throat_data Pc_drain_max Pc_interval n_p n_t
global sig_ow
sig_ow = 20e-3; % N/m

%% Calculating the pistone like threshold pressure for pores and throats
for ii = 1:n_t
    throat_data(ii,34) = throatThresholdPressPiston (throat_data(ii,18),...
        throat_data(ii,7),throat_data(ii,6),throat_data(ii,13),throat_data(ii,8:11));
end
for ii = 1:n_p
    pore_data(ii,28) = poreThresholdPressPiston (pore_data(ii,16), pore_data(ii,7),...
        pore_data(ii,6),pore_data(ii,12),pore_data(ii,8:11));
end
%% determining the capillary pressure level interval
[Pc_interval , Pc_drain_max , simTimes] = PcInterval;
A = pore_data(:,28);B = throat_data(:,34);
C = zeros(length(A) + length(B),1);
C(1:length(A)) = A;C(length(A) + 1:end) = B;
[Nelements,Xcenters] = hist(C,200);
D = Xcenters(Nelements > 10);

%%
% The backbone identifies which pores/throats are oil filled at the start of drainage.
% The first row is pore/throat index, followed by 1 for pores and 0 for thoats.
backbone = [];
fill = []; % A list for sequense of pore or throat fillings

t = 1;
Pc = 1;
while Pc < Pc_drain_max
    %% Cheking reservoir pore for potential of filling
    for ii = 1:n_p
        if pore_data(ii,3) == -1
            if pore_data(ii,28) <= Pc && pore_data(ii,5) == 0;
                pore_data(ii,5) = 1;
                fill(end + 1,1:2) = [ii , 1];
            end
        end
    end
    
    pore_list = [];
    throats_of_oil_invaded_pores = [];
    for i =1:n_p
        if  pore_data(i,5) == 1
            pore_list(i,1) = pore_data(i,1);
        end
    end
    if ~isempty(pore_list)
        pore_list(pore_list(:,1) == 0) = [];
        for i = 1:n_t
            if pore_data(throat_data(i,21),5) == 1 || pore_data(throat_data(i,22),5) == 1
                throat_list(i,1) = throat_data(i,1);
            end
        end
            throat_list(throat_list(:,1) == 0) = []; 

        while ~isempty(throat_list)
            [a,IX] = sort(throat_data(throat_list,34));
            throat_list = throat_list(IX);
            
            if throat_data(throat_list(1),34) <= Pc &&...
                    throat_data(throat_list(1),5) == 0 &&...
                    (pore_data(throat_data(throat_list(1),21),5) == 1 ||...
                    pore_data(throat_data(throat_list(1),22),5) == 1)

                throat_data(throat_list(1),5) = 1; % occupying the throat with oil
                fill(end + 1,1:2) = [throat_data(throat_list(1),1) , 0];
                
                if pore_data(throat_data(throat_list(1),21),5) ~= 1 &&...
                        pore_data(throat_data(throat_list(1),21),28) <= Pc

                    pore_data(throat_data(throat_list(1),21),5) = 1;
                    fill(end + 1,1:2) = [pore_data(throat_data(throat_list(1),21),1) , 1];
                    
                    pore_list(end + 1 ,1) = pore_data(throat_data(throat_list(1),21),1);
                    for i = 1:length(nonzeros(pore_data(pore_list(end),41:end)))
                        throat_list(end + 1 ,1) = nonzeros(pore_data(pore_list(end),40 + i));
                    end
                elseif pore_data(throat_data(throat_list(1),22),5) ~= 1 &&...
                        pore_data(throat_data(throat_list(1),22),28) <= Pc

                    pore_data(throat_data(throat_list(1),22),5) = 1;
                    fill(end + 1,1:2) = [pore_data(throat_data(throat_list(1),22),1) , 1];

                    pore_list(end + 1 ,1) = pore_data(throat_data(throat_list(1),22),1);
                    for i = 1:length(nonzeros(pore_data(pore_list(end),41:end)))
                      throat_list(end + 1 ,1) = nonzeros(pore_data(pore_list(end),40 + i));
                    end
                end
            end
            throat_list(throat_list(:) == throat_list(1)) = [];
        end
    end
    %% Updating element saturations and conductances
    for ii = 1:n_p
        if pore_data(ii,5) == 1
            % Piri
%             [pore_data(ii,27),pore_data(ii,26),pore_data(ii,22),pore_data(ii,21)] =...
%                 conduct_drain_4(Pc,pore_data(ii,16),pore_data(ii,8:11),...
%                 pore_data(ii,7),pore_data(ii,6),pore_data(ii,14),sig_ow,...
%                 pore_data(ii,12));
            % Valvatne
%             [pore_data(ii,27),pore_data(ii,26),pore_data(ii,22),...
%                 pore_data(ii,21)] = conduct_drain_3_1(Pc,pore_data(ii,16),...
%                 pore_data(ii,8:11),pore_data(ii,7),pore_data(ii,6),...
%                 pore_data(ii,14),sig_ow);
            %Patzek
            [pore_data(ii,27),pore_data(ii,26),pore_data(ii,22),pore_data(ii,21)] =...
                conduct_drain(Pc,pore_data(ii,16),pore_data(ii,8:11),...
                pore_data(ii,7),pore_data(ii,6),pore_data(ii,14),sig_ow,...
                pore_data(ii,12));
            
            % Fluid Volume Calculation
            [pore_data(ii,24),pore_data(ii,23)] = Cal_fluid_vol(pore_data(ii,22),...
                pore_data(ii,21),pore_data(ii,13));
        end      
    end
    for ii = 1:n_t
        if throat_data(ii,5) == 1
            % Piri
%             [throat_data(ii,33),throat_data(ii,32),throat_data(ii,28), throat_data(ii,27)] =...
%                 conduct_drain_4(Pc,throat_data(ii,18),throat_data(ii,8:11),throat_data(ii,7),...
%                 throat_data(ii,6),throat_data(ii,15),sig_ow,throat_data(ii,13));
            % Valvatne
%             [throat_data(ii,33),throat_data(ii,32),throat_data(ii,28), throat_data(ii,27)] =...
%                 conduct_drain_3_1(Pc,throat_data(ii,18),throat_data(ii,8:11),throat_data(ii,7),...
%                 throat_data(ii,6),throat_data(ii,15),sig_ow);
            %Patzek
            [throat_data(ii,33),throat_data(ii,32),throat_data(ii,28), throat_data(ii,27)] =...
                conduct_drain(Pc,throat_data(ii,18),throat_data(ii,8:11),throat_data(ii,7),...
                throat_data(ii,6),throat_data(ii,15),sig_ow,throat_data(ii,13));

            % Fluid Volume Calculation
            [throat_data(ii,30),throat_data(ii,29)] = Cal_fluid_vol(throat_data(ii,28),...
                throat_data(ii,27),throat_data(ii,14));
        end
    end    
    %% Preparing Pc , Sw & Kr data
    Pc_drain_curve(t,1) = Pc;
    
    % Relative Permeability Calculation
     [kr_oil(t,1),kr_water(t,1)] = k_rel(1,0);
    
     % Water Saturation Calculation
    water_vol = sum(pore_data(:,23)) + sum(throat_data(:,29));
    Sw_drain(t,1) = water_vol / totalPV;
    
    % Pc Step Calculation
    
    t = t + 1
    if t < length(D)
        Pc = D(t);
    else
        Pc = Pc + 5*Pc_interval;
    end
end
Pc_drain_curve = Pc_drain_curve*0.000145037738;


%% Ploting Pc & Kr
figure('name','Primary Dranage Cappilary Pressure & Relative Permeability Curves',...
    'units','normalized','outerposition',[0 0 1 1])
subplot(2,1,1);
grid on
plot(Sw_drain,Pc_drain_curve,'--r')
title('Drainage Cappilary Pressure Curves')
xlabel('Sw')
xlim([0 1.05])
ylabel('Pc (Pa)')
legend('Drainage Pc')

subplot(2,1,2);
plot(Sw_drain,kr_water,'-r',Sw_drain,kr_oil,'--b')
title('Drainage Relative Permeability Curves')
xlabel('Sw')
xlim([0 1.05])
ylabel('Reative Permeability')
ylim([0 1.05])
legend('Water Relative Permeability','Oil Relative Permeability','Location','West')

%%
save network.mat