clc; clear
fileName = 'Berea';

% Pressure difference
inletPressure = 1;
outletPressure = 0;

fluids = Fluids();
network = Network(fileName);
network.calculatePorosity();
fprintf('Porosity of the model is: %3.5f \n', network.Porosity)

network.pressureDistribution(1,0)
network.AbsolutePermeabilityCalculation();