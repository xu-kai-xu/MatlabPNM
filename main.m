clc; clear
fileName = 'Berea';

% Pressure difference
inletPressure = 1;
outletPressure = 0;
sig_ow = 0.02; 
water_viscosity = 0.001;
oil_viscosity = 0.001;
Pc = 3000;
network = Network(fileName);
fprintf('Porosity of the model is: %3.5f \n', network.calculatePorosity())

