function advancingAngleOfpores = poreAdvancingAngle (minporeAdvancingAngle ,...
    maxporeAdvancingAngle , delta , etha)
%% 
% Thsi function detremines pore advancing angle by the use of distribution 
% function called Weibull
%%
global n_p

x = rand(n_p,1);
advancingAngleOfpores = weibull(deg2rad(maxporeAdvancingAngle) ,...
    deg2rad(minporeAdvancingAngle) , delta , etha , x);