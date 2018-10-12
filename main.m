clc; clear
fileName = 'Berea';

% Pressure difference
inletPressure = 1;
outletPressure = 0;
network = Network(fileName);
fprintf('Porosity of the model is: %3.5f \n', network.calculatePorosity())

