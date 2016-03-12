function totalConductance = throatConductance(elementArea , ...
    elementShapeFactor , elementGeometry) 

global n_t water_viscosity
for ii = 1:n_t
    if elementGeometry(ii) == 1
        totalConductance(ii) = 0.5 * elementArea(ii)^2 * elementShapeFactor(ii) * water_viscosity;
    elseif elementGeometry(ii) == 2
        totalConductance(ii) = 3 * elementArea(ii)^2 * elementShapeFactor(ii) /water_viscosity / 5;
    elseif elementGeometry(ii) == 3
        totalConductance(ii) = 0.5623 * elementArea(ii)^2 * elementShapeFactor(ii) * water_viscosity;
    end
end