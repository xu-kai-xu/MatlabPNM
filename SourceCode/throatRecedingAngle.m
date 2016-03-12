function recedingAngleOfThroats = throatRecedingAngle (minThroatRecedingAngle ,...
    maxThroatRecedingAngle , delta , etha)
%% 
% Thsi function detremines throat receding angle by the use of distribution 
% function called Weibull
%%
global n_t

x = rand(n_t,1);
recedingAngleOfThroats = weibull(deg2rad(maxThroatRecedingAngle) ,...
    deg2rad(minThroatRecedingAngle) , delta , etha , x);