%% compare test and train data samples to validate decoding

model    = 'direction';
useRevco = false;

switch model
    case 'direction'
        baseDir='/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/residuals/data';
    case 'choice'
        switch useRevco
            case true
                baseDir='/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/residuals/data/allRevco';
            case false
                baseDir='/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/residuals/data';
        end
end

sList = dir(baseDir);
sList = sList(~ismember({sList.name},{'.','..','.DS_Store', 'oldBin', 'allRevco'})); % CLEAN UP

for iSession = 1:length(sList)
    
    % load full session
    session = load([baseDir filesep sList(iSession).name]);
    if isfield(session, 'session')
        session = session.session;
    end
    
    loTrialPredict  = cell2mat(arrayfun(@(x) x.isPredictionSame_trial , session.sLoo, 'UniformOutput', false));
    loTrialAvgAccuracy(iSession) = mean(loTrialPredict);

    allTrialPredict = cell2mat(arrayfun(@(x) x.isPredictionSame_all , session.sLoo, 'UniformOutput', false));
    loAllAvgAccuracy(iSession) = mean(allTrialPredict);
  
    trueAccuracy_all(iSession)   = session.sLoo(1).accuracyAll;
    trueAccuracy_train(iSession) = mean(cell2mat(arrayfun(@(x) x.accuracyTrain , session.sLoo, 'UniformOutput', false)));
end

%%

%figure
subplot(141); hold on
histogram(loTrialAvgAccuracy, 10)
xlabel('avg prediction accuracy on left-out trial')
ylabel('n sessions')

subplot(142); hold on
histogram(loAllAvgAccuracy, 10)
xlabel('avg prediction accuracy on all trials')
ylabel('n sessions')

subplot(143); hold on
plot(trueAccuracy_all - trueAccuracy_train, 'o'); hold on
plot([0 72], [0 0], 'k--')

subplot(144); hold on
plot(trueAccuracy_all, 'o')