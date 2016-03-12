% clc; clear 
load network.mat
global pore_data throat_data

%% Throat Data Construction
[pore_data , throat_data] = rawNetworkConstruction(numX , numY , numZ ,...
    minThroatR , maxThroatR , etaThroatR , deltaThroatR , ...
    minThroatL , maxThroatL , etaThroatL , deltaThroatL , ...
    minAspectRatio , maxAspectRatio , etaAspectRatio , deltaAspectRatio);
[throat_data(:,6) , throat_data(:,7)] = ...
    throatGeometry(squarePercThroat , circPercThroat , maxShapeFactThroat ,...
    minShapeFactThroat , etaShapeFactthroat , deltaShapeFactThroat);
[throat_data(:,8) , throat_data(:,9) , throat_data(:,10) , throat_data(:,11)] = ...
    throatAngles (throat_data(:,6));
[throat_data(:,15) , throat_data(:,16)] = ...
    throatAandV(throat_data(:,13) ,throat_data(:,6) , throat_data(:,14) , throat_data(:,7));
throat_data(:,18) = throatRecedingAngle (minRecedingAngle ,maxRecedingAngle , 0.1 , 1);
throat_data(:,19) = throatAdvancingAngle (minAdvancingAngle , maxAdvancingAngle , .4 , 1.5);
throat_data(:,31) = throatConductance(throat_data(:,15) , ...
    throat_data(:,6) , throat_data(:,7));
throat_data(:,17) = clayVolume;
throat_data(:,29) = throat_data(:,16); throat_data(:,29);throat_data(:,27) = throat_data(:,15);
throat_data(:,32) = throat_data(:,31);
throat_data(:,42) = nan; % holds cuslter Pc after imbibition
throat_data(:,4) = throat_data(:,5);
%% Pore Data Construction
[pore_data(:,6) , pore_data(:,7)] = poreGeometry(squarePercPore , circPercPore ,...
    maxShapeFactPore , minShapeFactPore , etaShapeFactPore , deltaShapeFactPore);
[pore_data(:,8) , pore_data(:,9) , pore_data(:,10) , pore_data(:,11)] =...
    poreAngles (pore_data(:,6));
[pore_data(:,14) , pore_data(:,15)] = ...
    poreAandV(pore_data(:,12) ,pore_data(:,6) , pore_data(:,13) , pore_data(:,7));
pore_data(:,16) = poreRecedingAngle (minRecedingAngle , maxRecedingAngle , 0.1 , 1);
pore_data(:,17) = poreAdvancingAngle (minAdvancingAngle , maxAdvancingAngle , 0.4 , 1.4);
pore_data(:,25) = ...
    poreConductance(pore_data(:,14) , pore_data(:,6) , pore_data(:,7));
pore_data(:,18) = clayVolume;
pore_data(:,21) = pore_data(:,14); pore_data(:,23) = pore_data(:,15);
pore_data(:,26) = pore_data(:,25);
pore_data(:,36) = nan; % hold cluster Pc after imbibitions 
pore_data(:,4) = pore_data(:,5);
%% Pore Volume Calculation
totalPV = totalPoreVol ;
%% Calculating the abslute permeability
abs_k_m2 = k_abs(1,0);
abs_k_d = abs_k_m2/9.869e-13;
fprintf('The permeability of the model is %2.4e m^2 which is %2.4f darcy\n',abs_k_m2,abs_k_d)
fprintf('Which is %3.2f milidarcy \n',abs_k_d *1000 )
fprintf('\n')
%% Visualozation Data
pore_data = poreLocation(pore_data);
modSize(1) = max(pore_data(:,37));
modSize(2) = max(pore_data(:,38));
modSize(3) = max(pore_data(:,39));
poreLoc = pore_data(:,37:39);
throatConn = throat_data(:,21:22);
%% Saving the constructed network
save network.mat