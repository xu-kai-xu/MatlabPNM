clc; clear
global oil_viscosity water_viscosity clayVolume

% Number of Nodes in different dimensions
numX = 3;
numY = 7;
numZ = 3;
% Throat radius data
minThroatR = 1e-06; maxThroatR = 20e-06; 
etaThroatR = 3; deltaThroatR = 0.2;
% Throat Length data
minThroatL = 20e-6; maxThroatL = 50e-6; 
etaThroatL = 3; deltaThroatL = 0.2;
% Aspect ratio data
minAspectRatio = 0.5; maxAspectRatio = 2; 
etaAspectRatio = .1; deltaAspectRatio = .5;
% Receding Angle data
minRecedingAngle = 0; maxRecedingAngle = 0;
% Advancing Angle data
minAdvancingAngle = 0; maxAdvancingAngle = 60;
% Percentage of each geometry
squarePercThroat = 15; circPercThroat = 85;
squarePercPore = 10; circPercPore = 85;
% Shape factor data
minShapeFactThroat = 0.01; maxShapeFactThroat = sqrt(3)/36;
etaShapeFactthroat = 1.8; deltaShapeFactThroat = 0.8;
minShapeFactPore = 0.01; maxShapeFactPore = sqrt(3)/36;
etaShapeFactPore = 1.8; deltaShapeFactPore = 0.6;
clayVolume = 0.2;
water_viscosity = 1e-3; %Pa_s = 1 cp ----(1e-3 = 1 cp)
oil_viscosity = 1e-3; %Pa_s = 1 cp ----(1e-3 = 1 cp)

save network.mat
NetworkConstruction
PrimaryDrainage
SecondaryImbibition