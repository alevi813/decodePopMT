
function [sessionStruct] = newDecodingTesters(model, condition, fitWindow, plotWindow)
% train a linear decoder on trial firing rates per neuron and the either
% the choice or the direction on each trial. Get weights for each neuron,
% reweight to get ensemble activity that best predicts the variable of your
% choice
% INPUTS:
% model     --- 'choice' or 'direction'
% condition --- entry {1} should be subj: 'nancy' or 'leo'
%           --- entry {2} is stimulus condition: 'flat', 'late', 'early'
% OUPUTS:
% sessionStruct with everything you should need...
% also automatically plots and saves some stuff right now. Come back later
% to put in full options for drawing + saving. Maybe add the ability to
% input specific session too?

%% get expts and set paths based on machine
experiments = getExperimentsAnd(condition);

comp = getComp;

if strcmp(comp, 'desktop')
    dataPath  = ['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/' condition{1}];
    nDataPath = ['/Users/Aaron/Dropbox/twagAnalysis4.1/Data/' condition{1} '/neurons'];
    
    if strcmp(model, 'choice')
        savePath   = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/choice/' condition{2}];
    else
        savePath   = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/direction/' condition{2}];
    end
else
    dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/' condition{1}];
    nDataPath = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/' condition{1} '/neurons'];
    
    if strcmp(model, 'choice')
        savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/sessionData/choice/' fitWindow 'Fit_' plotWindow 'Plot'];
    else
        savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/sessionData/direction/' fitWindow 'Fit_' plotWindow 'Plot'];
    end
end

nFiles = dir(nDataPath); %nFiles = nFiles(4:end);
nFiles=nFiles(~ismember({nFiles.name},{'.','..','.DS_Store'})); % clean up

%%

% get a list of all neuron file names
nFilesNames = arrayfun(@(x) x.name(1:9) , nFiles, 'UniformOutput', false);

for kEx = 1:numel(experiments)
    %for kEx = 2
    
    exname=experiments{kEx};
    %clear neurons
    % get neurons for the given session
    nNeurons = sum(strcmp(nFilesNames, exname));
    
    for iN = 1:nNeurons
        if iN < 10
            neurons(iN) = load([nDataPath filesep exname '_0' num2str(iN)]);
        else
            neurons(iN) = load([nDataPath filesep exname '_' num2str(iN)]);
        end
    end
    
    % load stim file
    stim    = getStim(exname, dataPath);
    S = behavior.sessionBehaviorSummary(stim);
    
    froCoh = S.Xc(S.froIx,:);
    allFroCoh(kEx, :) = froCoh(1,:);
    
    % stimulus timing variables, bin size, windows
    motionOnset = [stim.timing(:).motionon] + [stim.timing(:).plxstart];
    goTime      = [stim.timing(:).fpoff] + [stim.timing(:).plxstart];
    binSize = 0.01; % 10ms (in seconds)
    window  = [-.5 2]; % 500 ms before motion onset to ~1s after (in seconds)
    %window  = [-2.5 0]; % backwards windown when aligned to 'go'
    
    nBins = diff(window)/binSize;
    
    % spike rate for the population
    Rs = nan(stim.nTrials, nBins, nNeurons); %preallocate
    
    % loop over neurons in the session
    for kNeuron = 1:nNeurons
        
        spikeTimes  = neurons(kNeuron).spikeTimes;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% get spikes aligned to MOTION ONSET
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [motionOn.spcnt, motionOn.bins]  = pdsa.binSpTimes(spikeTimes, motionOnset, window, binSize);
        
        % smooth spike count with a boxcar filter
        sm = 5; % size of boxcar (number of bins)
        tmp = filter(boxcar(sm)/sm, 1, motionOn.spcnt');
        tmp = flipud(tmp);
        tmp = filter(boxcar(sm)/sm, 1, tmp);
        tmp = flipud(tmp);
        
        motionOn.spcnt = tmp';
        
        % conver to spike rate
        motionOn.sprate = motionOn.spcnt/binSize;
        
        goodTrials  = stim.goodtrial & ~any(isnan(motionOn.spcnt), 2);
        
        % get residuals (subtract off the effect of motion onset)
        motionOn.residuals                  = bsxfun(@minus, motionOn.sprate, nanmean(motionOn.sprate));
        motionOn.Rs(goodTrials, :, kNeuron) = motionOn.residuals(goodTrials,:);
        %         % using raw rates, not residuals
        %         motionOn.Rs(goodTrials, :, kNeuron) = motionOn.sprate(goodTrials,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% get spikes aligned to GO SIGNAL
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [fromGo.spcnt, fromGo.bins]  = pdsa.binSpTimes(spikeTimes, goTime, [-1 0], binSize);
        
        % smooth spike count with a boxcar filter
        sm = 5; % size of boxcar (number of bins)
        tmp_fromGo = filter(boxcar(sm)/sm, 1, fromGo.spcnt');
        tmp_fromGo = flipud(tmp_fromGo);
        tmp_fromGo = filter(boxcar(sm)/sm, 1, tmp_fromGo);
        tmp_fromGo = flipud(tmp_fromGo);
        
        fromGo.spcnt = tmp_fromGo';
        
        fromGo.sprate = fromGo.spcnt/binSize;
        
        %         % residuals
        %         fromGo.residuals                  = bsxfun(@minus, fromGo.sprate, nanmean(fromGo.sprate));
        %         fromGo.Rs(goodTrials, :, kNeuron) = fromGo.residuals(goodTrials,:);
        %         fromGo.spikes                     = fromGo.Rs(goodTrials,:,:);
        
        % using raw rates, not residuals
        fromGo.Rs(goodTrials, :, kNeuron) = fromGo.sprate(goodTrials,:);
        fromGo.spikes                    = fromGo.Rs(goodTrials,:,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end % neuron loop
    
    %% Now do population decoding approach
    
    % Rs is trialxbinsxneuron - sum across time bins, and then squeeze.
    % fitWindow determines whether you fit on stimulus (motionOn) or
    % nonstimulus time (fromGo)
    switch fitWindow
        case 'stimulus'
            % summed spikes over stimulus period +150ms
            R = squeeze(sum(motionOn.Rs(:,motionOn.bins > 0 & motionOn.bins < 1.2,:),2));
        case 'nonstimulus'
            % summed spikes over 500ms pre "go"
            R = squeeze(sum(fromGo.Rs(:,fromGo.bins > -.5 & fromGo.bins < 0,:),2));
    end
    
    Cho = sign(stim.targchosen - 1.5);
    Direc = sign(sum(sum(stim.pulses, 3),2));
    
    %goodTrials = ~any(isnan(R),2); %extraneous???
    
    R = R(goodTrials,:);
    Cho = Cho(goodTrials);
    Direc = Direc(goodTrials);
    dirprob = stim.dirprob(goodTrials);
    
    % which spikes are you projecting onto your weights? stimulus or
    % nonstimulus?
    switch plotWindow
        case 'stimulus'
            spikes = motionOn.Rs(goodTrials,:,:);
        case 'nonstimulus'
            spikes = fromGo.Rs(goodTrials,:,:);
    end
        
    if corr(Cho, Direc) < 0
        Cho = -Cho;
    end
    
    froIx = stim.trialId==stim.frozenTrialIds;
    froIx = froIx(goodTrials);
    
    % get weights using GLM fit    
    if strcmp(model, 'choice')
        % fit on choice, using *only frozen* trials
        % xval via leave-one-out
        w2 = glmfit( R(froIx, :), Cho(froIx) );
        [sLoo] = xval_LOO( R(froIx, :), Cho(froIx) );
        
        %         % fit on choice, using *all revco* trials
        %         w2 = glmfit( R(dirprob==0, :), Cho(dirprob==0) );
    else
        w2  = glmfit( R(~froIx, :), Direc(~froIx) ); % fit on direction (binary for now)
        [sLoo] = xval_LOO( R(~froIx, :), Direc(~froIx) );
        
    end % if choice or dir
    
    wAll = w2(2:end);
    %kAll = w2(1);
            
    if corr(R*wAll, Cho) < 0
        flipDV = true;
    else
        flipDV = false;
    end    
        
    %%%%%%%%%%%% get a decision variable for each trial %%%%%%%%%%%%%%
    %%% spikes var comes from switch statement above w/ plotWindow %%%
    [nTrials, nBins, nNeurons] = size(spikes);
    
    if strcmp(model, 'choice')
        dvAll   = nan(sum(froIx), nBins);
    else
        dvAll   = nan(nTrials, nBins);
    end
    
    % reweighted decision variable
    if strcmp(model, 'choice')
        froSpikes = spikes(froIx,:,:);
        for iTrial = 1:sum(froIx)
            R = squeeze(froSpikes(iTrial,:,:) ); %squeeze to get spikes per bin per neuron during a specific trial
            %dvAll(iTrial,:)   = R*wAll; % muliply spike matrix (nBin x nNeuron) by weight vector (nNeurons x 1)
            dvAll(iTrial,:)   = R* sLoo(iTrial).wTrain; % muliply spike matrix (nBin x nNeuron) by weight vector (nNeurons x 1)
            
            if flipDV
                dvAll(iTrial,:)   = -dvAll(iTrial,:);
            end
        end% ntrials
    else
        for iTrial = 1:nTrials
            R = squeeze(spikes(iTrial,:,:) ); %squeeze to get spikes per bin per neuron during a specific trial
            dvAll(iTrial,:)   = R*wAll; % muliply spike matrix (nBin x nNeuron) by weight vector (nNeurons x 1)
            %dvAll(iTrial,:)   = R* sLoo(iTrial).wTrain; % muliply spike matrix (nBin x nNeuron) by weight vector (nNeurons x 1)
            
            if flipDV
                dvAll(iTrial,:)   = -dvAll(iTrial,:);
            end
        end% ntrials
    end % choice v direction model
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % smooth the dv
    dvAll   = filter(ones(5,1)/5, 1, dvAll')';
    
    % calculate traditional cp
    % Cho is already aligned to pref direction of neuron above, but you
    % still need to binarize
    binaryCho = Cho==1;
    
    if strcmp(model, 'direction')
        [cp_m,cp_s] = choiceProbabilityCalculate(dvAll(froIx, :), binaryCho(froIx));
    else
        [cp_m,cp_s] = choiceProbabilityCalculate(dvAll, binaryCho(froIx));
    end
    
    
    
    %% package and save
    
    sessionStruct.dvAll     = dvAll;
    %     sessionStruct.dvTrain   = dvTrain;
    %     sessionStruct.dvTest     = dvTest;
    
    if strcmp(plotWindow, 'stimulus')
        sessionStruct.bins       = motionOn.bins;
    else 
        sessionStruct.bins       = fromGo.bins;
    end
    
    sessionStruct.spikes     = spikes;
    %sessionStruct.newdvAll   = [dvAll(binaryCho==1,:); (dvAll(binaryCho==0,:) .*-1)];
    sessionStruct.cho       = binaryCho;
    sessionStruct.direction = Direc;
    sessionStruct.froIx     = froIx;
    sessionStruct.cp_m      = cp_m;
    sessionStruct.cp_s      = cp_s;
    sessionStruct.stim      = stim;
    sessionStruct.behavior  = S;
    sessionStruct.wAll      = wAll;
    sessionStruct.sLoo      = sLoo;
    sessionStruct.flipDV    = flipDV;
    sessionStruct.nNeurons  = nNeurons;
    sessionStruct.fitWindow = fitWindow;
    sessionStruct.plotWindow = plotWindow;
    
    % save the session data strcut
    save([savePath filesep exname], '-v7.3', '-struct', 'sessionStruct')
 
    clearvars -except model condition weightOrigin savePath experiments comp nFiles nFilesNames nDataPath dataPath fitWindow plotWindow
end % experiment loop


