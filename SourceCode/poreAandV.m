function [poreArea , poreVolume] = poreAandV(poreRadius ,...
    poreShapeFactor , poreLength , poreGeometry)
%%
% This function calculates the cross section area and volume of each pores
% (pore or pores)
%%
global n_p clayVolume
for ii = 1:n_p
    if poreGeometry(ii) == 1
        poreArea(ii) = (1+clayVolume)*( pi * poreRadius(ii) ^ 2);
        poreVolume(ii) = poreArea(ii) * poreLength(ii); 
    elseif poreGeometry(ii) == 2
        poreArea(ii) = (1+clayVolume)*( poreRadius(ii) ^ 2 / 4 / poreShapeFactor(ii));
        poreVolume(ii) = poreLength(ii) * poreArea(ii) ;
    elseif poreGeometry(ii) == 3
        poreArea(ii) = (1+clayVolume)*( 4 * poreRadius(ii) ^ 2);
        poreVolume(ii) = poreLength(ii) * poreArea(ii) ;    
    end
end

