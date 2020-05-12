% % direction
% session = load('/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/residuals/data/l20190806.mat');
% % choice
% %session = load('/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/residuals/data/l20190806.mat');


condition = 'early';
model     = 'direction';
experiments = getExperimentsAnd(condition);
dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model '/residuals/data'];

combo.sumcoh       = [];
combo.pulsecoh     = [];
combo.isRevco      = [];
combo.neuralChoice = [];

for iExp = 1:length(experiments)
    session = load([dataPath filesep experiments{iExp} '.mat']);
    
    if session.nNeurons > 1
    spikeLeft  = session.dvAll(session.direction==-1,:);
    spikeRight = session.dvAll(session.direction==1,:);
    
    spikeLeft  = sum(spikeLeft(:, 50:170),2);
    spikeRight = sum(spikeRight(:, 50:170),2);
    
    spikesAll = sum(session.dvAll(:, 50:170), 2);
    
    pulsecoh = sum(session.stim.pulses, 3);
    sumcoh   = mean(pulsecoh, 2) ./19;
    isRevco  = session.stim.dirprob == 0;
    
    if isfield(session.stim, 'validTrials')
        pulsecoh = pulsecoh(session.stim.validTrials, :);
        sumcoh   = sumcoh(session.stim.validTrials);
        isRevco  = isRevco(session.stim.validTrials);
    else
        pulsecoh = pulsecoh(session.stim.goodtrial, :);
        sumcoh   = sumcoh(session.stim.goodtrial);
        isRevco  = isRevco(session.stim.goodtrial);
    end
    
    neuralChoice = spikesAll>0;
    %neuralChoice = spikesAll> -0.1355;
    
    %%
    figure
    subplot(121)
    pmf  = pmfTools(sumcoh, neuralChoice, 'nBins', 18);
    fit(pmf);
    plot(pmf);
    
    subplot(122)
    ppk = ppkTools(pulsecoh, neuralChoice, 'ridge', true);
    plot(ppk, 'plotFit', false);
    
    %%
    combo.sumcoh       = [combo.sumcoh; sumcoh];
    combo.pulsecoh     = [combo.pulsecoh; pulsecoh];
    combo.isRevco      = [combo.isRevco; isRevco];
    combo.neuralChoice = [combo.neuralChoice; neuralChoice];
    end %nNeuron check
end % session loop

%% plot pmf and ppk using data from all sessions
figure(222)
subplot(325)
pmf  = pmfTools(combo.sumcoh, combo.neuralChoice, 'nBins', 18);
fit(pmf);
plot(pmf);
title(condition, 'fontsize', 12)


subplot(326)
ppk = ppkTools(combo.pulsecoh, combo.neuralChoice, 'ridge', true);
plot(ppk, 'plotFit', false); hold on

% ppk_rc = ppkTools(combo.pulsecoh(combo.isRevco==1,:), combo.neuralChoice(combo.isRevco==1), 'ridge', true);
% plot(ppk_rc, 'plotFit', false);


%supertitle([ condition ' - all sessions'], 12)