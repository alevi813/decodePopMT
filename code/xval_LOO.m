function [xvalStruct] = xval_LOO(x, y)
% Cross validate via 'leave one out.' Made for twag decoding.
% x is trialx1 firing rate, y is binary var (choice or dir), dataSize is
% the proportion of data you want to leave out.

nTrials = size(y, 1);

% fit the weights, get reweighted dv
wAll    = glmfit(x, y); wAll = wAll(2:end);
dvAll   = x*wAll;

% get a binarized choice/direction var
y2 = y;
y2(y2==-1) = 0;

% make predictions based on dv, measure the overall accuracy
predictAll  = dvAll>0;
accuracyAll = sum(predictAll==y2) / length(y2);

% Leave out a trial, then fit weights.
% Just do it in order, cause that's easiest.
for kFold = 1:nTrials
    % set up indices
    allT     = 1:nTrials;
    trainIx  = find(allT~=kFold);
    testIx   = kFold;
    
    % get new x&y vals
    xTrain   = x(trainIx, :);
    yTrain   = y(trainIx,:);
    
    xTest   = x(testIx, :);
    yTest   = y(testIx,:);
    
    % find weights, get reweighted dv
    tmpW  = glmfit(xTrain, yTrain); 
    xvalStruct(kFold).wTrain = tmpW(2:end);
    dvTrain = x* xvalStruct(kFold).wTrain;
    % make predictions from dv
    xvalStruct(kFold).predictTrain = dvTrain >0;
    
    % measure agreement between predictions of left-out trial and all
    % trials
    xvalStruct(kFold).isPredictionSame_all    = sum(xvalStruct(kFold).predictTrain==predictAll) / length(predictAll);
    xvalStruct(kFold).isPredictionSame_trial  = xvalStruct(kFold).predictTrain(kFold) == predictAll(kFold);
    
    % record the overall accruacy
    xvalStruct(kFold).accuracyTrain = sum(xvalStruct(kFold).predictTrain==y2)/ length(y2);
    
    % rename some stuff from early in each entry of the struct array
    xvalStruct(kFold).predictAll  = predictAll;
    xvalStruct(kFold).accuracyAll = accuracyAll;

end
