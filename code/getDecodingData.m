function [dS] = getDecodingData(condition)
% Load decoding data from a particular condition
%   condition: all, flat, late, or early

if nargin < 1
    condition = 'all';
end

if ~strcmp(condition, 'all')
    if strcmp(condition(2), 'all')
        condList = {'flat', 'late', 'early'};
        for iCond = 1:3
            dS.(condition{1}).(condList{iCond}).dir = decodingConditionAverage('direction', {condition(1), condList(iCond)});
            dS.(condition{1}).(condList{iCond}).cho = decodingConditionAverage('choice', {condition(1), condList(iCond)});
        end
    else
        dS.dir = decodingConditionAverage('direction', condition);
        dS.cho = decodingConditionAverage('choice', condition);
    end    
else
    if strcmp(condition, 'all')
        condList = {'flat', 'late', 'early'};
        for iCond = 1:3
            dS.(condList{iCond}).dir = decodingConditionAverage('direction', condList(iCond));
            dS.(condList{iCond}).cho = decodingConditionAverage('choice', condList(iCond));
        end
    else
        dS.dir = decodingConditionAverage('direction', condition);
        dS.cho = decodingConditionAverage('choice', condition);
    end
end

end

