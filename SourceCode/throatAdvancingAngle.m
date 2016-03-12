function advancingAngleOfThroats = throatAdvancingAngle (minthroatAdvancingAngle ,...
    maxthroatAdvancingAngle , delta , etha)
%% 
% Thsi function detremines throat advancing angle by the use of distribution 
% function called Weibull
%%
global n_t

x = rand(n_t,1);
advancingAngleOfThroats = weibull(deg2rad(maxthroatAdvancingAngle) ,...
    deg2rad(minthroatAdvancingAngle) , delta , etha , x);