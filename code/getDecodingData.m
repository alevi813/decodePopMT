function [dS] = getDecodingData(condition, fitWindow, plotWindow)
% Load decoding data from a particular condition
%   condition: all, flat, late, or early

if ~strcmp(condition, 'all')
    if strcmp(condition(2), 'all')
        condList = {'flat', 'late', 'early'};
        for iCond = 1:3
            dS.(condition{1}).(condList{iCond}).dir = decodingConditionAverage('direction', {condition(1), condList(iCond)}, fitWindow, plotWindow);
            dS.(condition{1}).(condList{iCond}).cho = decodingConditionAverage('choice', {condition(1), condList(iCond)}, fitWindow, plotWindow);
        end
    else
        dS.dir = decodingConditionAverage('direction', condition, fitWindow, plotWindow);
        dS.cho = decodingConditionAverage('choice', condition, fitWindow, plotWindow);
    end
else
    if strcmp(condition, 'all')
        condList = {'flat', 'late', 'early'};
        for iCond = 1:3
            dS.(condList{iCond}).dir = decodingConditionAverage('direction', condList(iCond), fitWindow, plotWindow);
            dS.(condList{iCond}).cho = decodingConditionAverage('choice', condList(iCond), fitWindow, plotWindow);
        end
    else
        dS.dir = decodingConditionAverage('direction', condition, fitWindow, plotWindow);
        dS.cho = decodingConditionAverage('choice', condition, fitWindow, plotWindow);
    end
end

end

