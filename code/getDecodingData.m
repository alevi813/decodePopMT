function [dS] = getDecodingData(condition, backward)
% Load decoding data from a particular condition
%   condition: all, flat, late, or early

if nargin <2
    backward = false;
    if nargin < 1
        condition = 'all';
    end
end

if ~strcmp(condition, 'all')
    if strcmp(condition(2), 'all')
        condList = {'flat', 'late', 'early'};
        for iCond = 1:3
            dS.(condition{1}).(condList{iCond}).dir = decodingConditionAverage('direction', {condition(1), condList(iCond)}, backward);
            dS.(condition{1}).(condList{iCond}).cho = decodingConditionAverage('choice', {condition(1), condList(iCond)}, backward);
        end
    else
        dS.dir = decodingConditionAverage('direction', condition, backward);
        dS.cho = decodingConditionAverage('choice', condition, backward);
    end
else
    if strcmp(condition, 'all')
        condList = {'flat', 'late', 'early'};
        for iCond = 1:3
            dS.(condList{iCond}).dir = decodingConditionAverage('direction', condList(iCond), backward);
            dS.(condList{iCond}).cho = decodingConditionAverage('choice', condList(iCond), backward);
        end
    else
        dS.dir = decodingConditionAverage('direction', condition, backward);
        dS.cho = decodingConditionAverage('choice', condition, backward);
    end
end

end

