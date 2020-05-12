%% example neurons
comp = getComp;

condition = 'flat';

% first plot scatter plots of decoding weights
dir = decodingConditionAverage('direction', condition);
cho = decodingConditionAverage('choice', condition);

% % try just a session instead of full dataset
% exname = 'l20191120';
%
% if strcmp(comp, 'laptop')
%     cho = load(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/data/' exname '.mat']);
%     dir = load(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/data/' exname '.mat']);
% else
%     cho = load(['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/choice/data/' exname '.mat']);
%     dir = load(['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/direction/data/' exname '.mat']);
% end

figure; hold on
%scatter(dir.wAll, cho.wAll)
plot(dir.w, cho.w, 'o') %, 'color', [255/255 179/255 0/255])
xlim([-2e-3 2e-3])
ylim([-3.5e-3 3.5e-3])
xlabel('wDir', 'FontSize', 14)
ylabel('wCho', 'FontSize', 14)

plot([0 0], [-3.5e-3 3.5e-3], 'k--')
plot([-2e-3 2e-3], [0 0], 'k--')
axis square

%% ehhhhh just load all neurons for now.

% set up paths
if strcmp(comp, 'laptop')
    dataPath{1} = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy';
    dataPath{2} = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/leo';
else
    dataPath{1} = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/nancy';    
    dataPath{2} = '/Users/Aaron/Dropbox/twagAnalysis4.1/Data/leo';
end

% load nancy neurons
%dataPath{1} = getpref('mtlipglm', 'dataPath');
fitDir    = 'fits_cho1'; % main_fits / fits_notargs / fits_pulsecovars / fits_cho1
S_nancy = loadUpFits(dataPath{1}, 'MT', condition, fitDir); % rate predictions under choice 1 kernel (right)

% load leo neurons
fitDir    = 'fits_leo'; % main_fits / fits_notargs / fits_pulsecovars / fits_cho1
S_leo = loadUpFits(dataPath{2}, 'MT', condition, fitDir); % rate predictions under choice 1 kernel (right)

% combine
S = [S_nancy S_leo];

% % take neurons only from relevant session, if that's what you want
% neuronExLabels = arrayfun(@(x) x.exname , S, 'UniformOutput', false);
% exIx = contains(neuronExLabels, exname);
% S = S(exIx);


%% quick diagnostic... compare d' and wDir
dp = arrayfun(@(x) x.model(1).dprime, S);

skipthis = true;
if ~skipthis
    figure
    subplot(1,2,1)
    scatter(dp, dir.w)
    %scatter(dp, dir.wAll)
    lsline
    xlabel('d''', 'fontsize', 14)
    ylabel('wDir', 'fontsize', 14)
    rval = corrcoef(dp, dir.w);
    %rval = corrcoef(dp, dir.wAll);
    %text(1, -.6e-3, ['r = ' num2str(rval(1,2))])
    title(['r = ' num2str(rval(1,2))])
    axis square
    %% cp vs. wCho
    
    cp = grandCP_test(S, dataPath);
    
    subplot(1,2,2)
    scatter(abs(cp-.5), cho.w)
    %scatter(abs(cp-.5), cho.wAll)
    lsline
    xlabel('CP', 'fontsize', 14)
    ylabel('wCho', 'fontsize', 14)
    rval = corrcoef(abs(cp-.5), cho.w);
    %rval = corrcoef(abs(cp-.5), cho.wAll);
    %text(1, -.6e-3, ['r = ' num2str(rval(1,2))])
    title(['r = ' num2str(rval(1,2))])
    axis square
end
%% single out some example neurons
%  some example cells from the flat condition:
%  +,+ cell 5
%  0,+ cell 155
%  +,0 cell 43

cellId = 5;

%for cellId = 1:length(dir.wAll)
thisNeuron = S(cellId);
% figure
% subplot(1,2,1);
% responseByDist(thisNeuron, dataPath);
% subplot(1,2,2);
% plotCPt(thisNeuron, dataPath);

%%
figure
% get dirprob (stimDistNum)

if strcmp(thisNeuron.exname(1), 'n')
    stim = load([dataPath{1} filesep 'stim' filesep thisNeuron.exname '_stim.mat']);
else
    stim = load([dataPath{2} filesep 'stim' filesep thisNeuron.exname '_stim.mat']);    
end

froIx = stim.trialId==stim.frozenTrialIds;
if isfield(stim, 'validTrials')
    froIx   = froIx(stim.validTrials);
    dirprob = stim.dirprob(stim.validTrials);
    targcorrect = stim.targcorrect(stim.validTrials);
    targchosen = stim.targchosen(stim.validTrials);
else
    froIx   = froIx(stim.goodtrial);    
    dirprob = stim.dirprob(stim.goodtrial);
    targcorrect = stim.targcorrect(stim.goodtrial);
    targchosen = stim.targchosen(stim.goodtrial);
end

%trialRates = reshape( cell2mat(arrayfun(@(x) x.model(1).trialRates/max(x.model(1).trialRates(:)), thisNeuron, 'UniformOutput', false)) , [size(thisNeuron.model(1).trialRates) 1]);

% rename trial rates
trialRates = S(cellId).model(1).trialRates;

% fix for sessions where I hit record after starting pds...
if length(trialRates) ~= length(dirprob)
    lDiff = length(dirprob)-length(trialRates);
    dirprob     = dirprob(lDiff+1:end);
    targcorrect = targcorrect(lDiff+1:end);
    targchosen  = targchosen(lDiff+1:end);
end

% direction psth
subplot(221); hold on
psthTime = S(1).model(1).psthTime;
% plot(psthTime , nanmean(trialRates(targcorrect==1,:)))
% plot(psthTime , nanmean(trialRates(targcorrect==2,:)))

dirRates1.m = nanmean(trialRates(targcorrect==1,:));
dirRates2.m = nanmean(trialRates(targcorrect==2,:));

dirRates1.se = nanstd(trialRates(targcorrect==1,:))  / sqrt(size(trialRates, 1));
dirRates2.se = nanstd(trialRates(targcorrect==2,:))  / sqrt(size(trialRates, 1));

% normalize?
dirRates1.mNorm = dirRates1.m / max(dirRates1.m);
dirRates2.mNorm = dirRates2.m / max(dirRates1.m);

boundedline(psthTime, smooth(dirRates1.m, 5), smooth(dirRates1.se, 5), 'cmap', [0 0 .6], 'alpha');
boundedline(psthTime, smooth(dirRates2.m, 5), smooth(dirRates2.se, 5), 'cmap', [.6 0 0], 'alpha');
%boundedline(psthTime, smooth(dirRates1.mNorm, 5), smooth(dirRates1.se, 5), 'cmap', [227/255 148/255 0], 'alpha');
%boundedline(psthTime, smooth(dirRates2.mNorm, 5), smooth(dirRates2.se, 5), 'cmap', [255/255 211/255 89/255], 'alpha');
title('direction PSTH -- all trials')
xlim([-500 1500])
axis square

% choice psth
cho1ix = froIx & targchosen==1;
cho2ix = froIx & targchosen==2;

choRates1.m = nanmean(trialRates(cho1ix,:));
choRates2.m = nanmean(trialRates(cho2ix,:));

choRates1.se = nanstd(trialRates(cho1ix,:)) / sqrt(size(trialRates(cho1ix,:), 1));
choRates2.se = nanstd(trialRates(cho2ix,:)) / sqrt(size(trialRates(cho1ix,:), 1));

subplot(222); hold on
% boundedline(psthTime, choRates1.m, choRates1.se, 'cmap', [0 .7 .7], 'alpha');
% boundedline(psthTime, choRates2.m, choRates2.se, 'cmap', [0 .4 .4], 'alpha');
boundedline(psthTime, smooth(choRates1.m, 5), smooth(choRates1.se, 5), 'cmap', [0 0 .6], 'alpha');
boundedline(psthTime, smooth(choRates2.m, 5), smooth(choRates2.se, 5), 'cmap', [.6 0 0], 'alpha');
title('choice PSTH -- frozen trials')
xlim([-500 1500])
axis square

% neurometric
[xAll, allAuc, aucse] = nmfByNeuron(thisNeuron, dataPath);
[xfit, fitvals]       = logFit(xAll, allAuc);

subplot(223)
errorbar(xAll, allAuc, aucse, 'o'); hold on
plot(xfit, fitvals)
axis square

subplot(224);
plotCPt(thisNeuron, dataPath);
axis square
set(gcf, 'color', 'white')

%end
% things to add?
%   - cp vs choW scatter
%   - choice sorted PSTH (raw)
%   - neurometric (raw)