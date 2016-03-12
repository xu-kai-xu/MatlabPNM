function [poreShapeFactor , poreGeometry] = poreGeometry(squarePorePercentage ,...
    tiangularPorePercentage , maxShapeFactorOfTriangularPores ,...
    minShapeFactorOfTriangularPores , delta , etha)
%%
% This function detrmines the percentage of each pore geometry (square, 
% triangular and circular geometries) and shape factor for each element.

global n_p
%% Triangular Section
x = rand(n_p,1);

if maxShapeFactorOfTriangularPores > sqrt(3)/36
    maxShapeFactorOfTriangularPores = sqrt(3)/36;
end

tiangularPorePercentage = tiangularPorePercentage/100;

tri_pore = x > (1 - tiangularPorePercentage);
poreGeometry = zeros(n_p,1);
poreShapeFactor = zeros(n_p,1);
poreGeometry(tri_pore) = 2; % 2 stands for triangular geometry

poreShapeFactor(tri_pore) = weibull(maxShapeFactorOfTriangularPores ,...
    minShapeFactorOfTriangularPores , delta , etha , rand(length(nonzeros(tri_pore)),1)); 

fprintf('%2.2f Percentage of the pores are Triangular.\n' ,size(nonzeros(tri_pore),1)/n_p*100);

%% Square Section
squarePorePercentage = squarePorePercentage / 100;
squ_pore = x <= (1 - tiangularPorePercentage)  & x > (1 - tiangularPorePercentage - squarePorePercentage);
poreGeometry(squ_pore) = 3; % 3 stands for sqular geometry
poreShapeFactor(squ_pore) = 0.0625;
fprintf('%2.2f Percentage of the pores are Square.\n' ,size(nonzeros(squ_pore),1)/n_p*100);

%% Circular Section
fprintf('%2.2f Percentage of the pores are Circular.\n' ,...
    (n_p - size(nonzeros(squ_pore),1)-size(nonzeros(tri_pore),1))/n_p*100);
for ii=1:length(poreGeometry)
    if poreGeometry(ii) == 0
        poreGeometry(ii) = 1; % 1 stands for circular geometry
        poreShapeFactor(ii) = 1/(4*pi);
    end
end
fprintf('\n')