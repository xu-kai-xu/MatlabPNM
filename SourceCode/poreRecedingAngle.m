function recedingAngleOfpores = poreRecedingAngle (minporeRecedingAngle ,...
    maxporeRecedingAngle , delta , etha)
%% 
% Thsi function detremines pore receding angle by the use of distribution 
% function called Weibull
%%
global n_p

x = rand(n_p,1);
recedingAngleOfpores = weibull(deg2rad(maxporeRecedingAngle) ,...
    deg2rad(minporeRecedingAngle) , delta , etha , x);