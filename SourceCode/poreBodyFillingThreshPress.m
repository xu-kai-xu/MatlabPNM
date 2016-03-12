function [poreBodyFillThreshPress , layerColapseThreshPress] =...
    poreBodyFillingThreshPress(poreIndex,poreLayer,advancingingAngle,...
    geometry,shapeFactor,inscribedR,halfAngles,rec_angle)

global pore_data throat_data sig_ow
W = [0;0.72;0.45;1.2;1.5;5];

attachedThroats = pore_data(poreIndex,41:end);
oilFilledAttachedThroats = zeros(length(nonzeros(attachedThroats)),1);
for i = 1:length(nonzeros(attachedThroats))
    if throat_data(attachedThroats(i),5) == 1
        oilFilledAttachedThroats(i,1) = attachedThroats(i);
    end
end
oilFilledAttachedThroats = nonzeros(oilFilledAttachedThroats);
z = length(oilFilledAttachedThroats); % number of oil filled attached throats

if poreLayer == -2
    [poreBodyFillThreshPress,layerColapseThreshPress] =...  
        pistonLikeThreshPressImb(geometry,halfAngles,rec_angle,...
        advancingingAngle,inscribedR,shapeFactor,sig_ow);
else
    if z == 0
        poreBodyFillThreshPress = nan;
        layerColapseThreshPress = nan;

    elseif z == 1
        [poreBodyFillThreshPress,layerColapseThreshPress] =...
            pistonLikeThreshPressImb(geometry,halfAngles,rec_angle,...
            advancingingAngle,inscribedR,shapeFactor,sig_ow);
    else
            if z > 5
                w = W(6);
            else
                w = W(z);
            end
            nominator = 0;
            denominator = 0;
            sumOfThroatRadius = 0;
            for ii = 1:z
                randNumber = rand;
                sumOfThroatRadius = sumOfThroatRadius + throat_data(oilFilledAttachedThroats(z),13);
                nominator = nominator + randNumber * sumOfThroatRadius;
                denominator = denominator + randNumber;
            end
            R_ave = (pore_data(poreIndex,12) + w * nominator / denominator)*cos(advancingingAngle);
            poreBodyFillThreshPress = 2*sig_ow/R_ave;

            layerColapseThreshPress = nan;
    end
end



        