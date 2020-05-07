% test validation

model       = 'direction';
condition   = 'early';
experiments = getExperimentsAnd({condition});

dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model '/data'];
figPath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model filesep condition];

doPlot = true;

for iExp = 1:length(experiments)
    % load session
    session = load([dataPath filesep experiments{iExp} '.mat']);
    
    
    [nTrials, nBins, nNeurons] = size(session.spikes);
    sessionNeurons(iExp) = nNeurons;
    
    testSpikes = session.spikes(session.testIx, :, :);
    testDir    = session.direction(session.testIx);
    noDirTest  = testDir==0;
    testDir    = testDir(~noDirTest)==1;
    testSpikes = testSpikes(~noDirTest, :, :);
    
    for iT = 1:nBins
        R = squeeze(testSpikes(:,iT,:) );
        dvVal(:,iT) = R*session.wTrain ;
        
        if session.flipDV
            dvVal(:,iT) = -dvVal(:,iT);
        end
    end
    
    % smooth the dv
    dvVal = filter(ones(5,1)/5, 1, dvVal')';
    
    auc(iExp) = roc(mean(dvVal(:,31:102), 2), testDir);
    
    if auc(iExp).AUC < .5
        auc(iExp).AUC = 1 - auc(iExp).AUC;
    end
    
    clear testSpikes testDir dvVal
end %exp

testAccuracy.m  = cell2mat(arrayfun(@(x) x.AUC, auc, 'UniformOutput', false));
testAccuracy.se = cell2mat(arrayfun(@(x) x.serror, auc, 'UniformOutput', false));

