function combo = decodingConditionAverage(model, condition, fitWindow, plotWindow)
% decodingConditionAverage
% INPUTS:
% - model --- 'choice' or 'direction'
% - condition --- can be {2x1} with {'monkey', 'condition'} or {1x1} with
% just condition for both monkeys
% OUTPUTS:
% - combo --- struct with aggregate vars concatenated from various
% individual session decoding files.


experiments = getExperimentsAnd(condition);

comp = getComp;

if strcmp(comp, 'desktop')
    if strcmp(model, 'direction')
        dataPath  = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/' model '/data'];
        %dataPath  = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/' model '/data/rawRates'];
    else
        %        dataPath  = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/' model '/data/allRevco'];
        %        dataPath  = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/' model '/data/rawRates'];
        dataPath  = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/' model '/data'];
    end
    figPath   = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/' model filesep condition];
else
    %     dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model '/rawRates/data'];
    %     figPath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model filesep 'rawRates/' condition];
    
%        dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model '/residuals/data'];
%        figPath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model filesep 'residuals/' condition];
    
    %     dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model '/xval_rawRates/data'];
    %     figPath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model filesep 'test_xval/' condition];
    
    
    %dataPath = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/sessionData/no_xval/' model filesep fitWindow 'Fit_' plotWindow 'Plot'];
    dataPath = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/sessionData/xvalMean/' model filesep fitWindow 'Fit_' plotWindow 'Plot'];
    %dataPath = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/sessionData/' model filesep fitWindow 'Fit_' plotWindow 'Plot'];
end


% combin the vars you want across sessions
combo.dvAll   = [];
combo.dvTrain = [];
combo.dvTest  = [];
combo.cho = [];
combo.pta = {};
combo.ppk = [];
combo.cpt = [];
combo.dir = [];
combo.spikes  = [];
combo.froIx   = [];
combo.dirprob = [];
combo.sumcoh  = [];
combo.targRchosen = [];
combo.nNeurons = [];
combo.w = [];
combo.slooW = [];
combo.timeToGo = [];

for iExp = 1:length(experiments)
    
    session = load([dataPath filesep experiments{iExp} '.mat']);
    
    if session.nNeurons > 1     % else do nothing
        pulses = sum(session.stim.pulses, 3);
        sumCoh = mean(pulses,2);
        
        if isfield(session.stim, 'validTrials')
            combo.dirprob = [combo.dirprob; session.stim.dirprob(session.stim.validTrials)];
        else
            combo.dirprob = [combo.dirprob; session.stim.dirprob(session.stim.goodtrial)];
        end
        
        tmpSpikes      = mean(session.spikes, 3);
        combo.spikes   = [combo.spikes; tmpSpikes];
        bins  = session.bins;
        
        slooWeights = cell2mat(arrayfun(@(x) x.wTrain , session.sLoo, 'UniformOutput', false));
        slooWeights = mean(slooWeights, 2);
        
        combo.sumcoh   = [combo.sumcoh; sumCoh];
        combo.froIx    = [combo.froIx; session.froIx];
        combo.cho      = [combo.cho; session.cho];
        combo.dvAll    = [combo.dvAll; session.dvAll];
        %        combo.dvTrain  = [combo.dvTrain; session.dvTrain];
        %        combo.dvTest   = [combo.dvTest; session.dvTest];
        %combo.pta      = [combo.pta; session.pta];
        combo.ppk      = [combo.ppk; (session.behavior.pk.mle(1:7) ./ norm(session.behavior.pk.mle(1:7)))'];
        combo.cpt      = [combo.cpt; session.cp_m'];
        combo.dir      = [combo.dir; session.direction];
        combo.nNeurons = [combo.nNeurons; session.nNeurons];
        combo.w        = [combo.w; session.wAll];
        combo.slooW    = [combo.slooW; slooWeights];
    end
    
    clear session
end % iExp loop

combo.bins = bins;
combo.pta  = cat(3, combo.pta{:});
combo.mPTA = mean(combo.pta, 3);

% newdv       = [combo.dvAll(combo.cho==1,:); (combo.dvAll(combo.cho==0,:) .*-1)];
% combo.newdv = newdv;

%normdv = mean(newdv) / max(mean(newdv));


% plot(bins, normdv); hold on

