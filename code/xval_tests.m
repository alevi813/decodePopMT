
function [sessionStruct] = xval_tests(model, condition)
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
        %savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/rawRates/' condition{2}];
        %savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/residuals/' condition{2}];
        savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/test_xval/' condition{2}];
    else
        %savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/rawRates/' condition{2}];
        %savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/residuals/' condition{2}];
        savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/test_xval/' condition{2}];
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
    
    %% Loop over neurons and plot their reponses
    figure(1); clf
    % plot pmf and kernel
    subplot(3,2,1)
    S.pmf.fit();
    S.pmf.plot();
    axis tight
    subplot(3,2,2)
    plot(S.pk.mle(1:7), 'k-o'); hold on
    plot(S.pk.mle(1:7)+S.pk.sd(1:7), '--k');
    plot(S.pk.mle(1:7)-S.pk.sd(1:7), '--k');
    xlabel('Pulse')
    ylabel('PPK Weight')
    axis tight
    
    %binSize = 1/60; % frame rate (in seconds)
    binSize = 0.01; % 10ms (in seconds)
    window  = [-.5 2]; % 500 ms before motion onset to ~1s after (in seconds)
    %window  = [-2.5 0]; % backwards windown when aligned to 'go'
    
    nBins = diff(window)/binSize;
    
    % spike rate for the population
    Rs = nan(stim.nTrials, nBins, nNeurons); %preallocate
    
    % loop over neurons in the session
    for kNeuron = 1:nNeurons
        
        
        spikeTimes  = neurons(kNeuron).spikeTimes;
        motionOnset = [stim.timing(:).motionon] + [stim.timing(:).plxstart];
        goTime      = [stim.timing(:).fpoff] + [stim.timing(:).plxstart];
        
        % get spike count aligned to motion onset
        [spcnt, bins]  = pdsa.binSpTimes(spikeTimes, motionOnset, window, binSize);
        %[spcnt, bins]  = pdsa.binSpTimes(spikeTimes, goTime, window, binSize);
        
        % smooth spike count with a boxcar filter
        sm = 5; % size of boxcar (number of bins)
        tmp = filter(boxcar(sm)/sm, 1, spcnt');
        tmp = flipud(tmp);
        tmp = filter(boxcar(sm)/sm, 1, tmp);
        tmp = flipud(tmp);
        
        spcnt = tmp';
        
        % conver to spike rate
        sprate = spcnt/binSize;
        
        goodTrials  = stim.goodtrial & ~any(isnan(spcnt), 2);
        
        
        % --- plot choice-sorted PSTH
        
        subplot(3,nNeurons,nNeurons + kNeuron)
        
        set(gcf, 'DefaultAxesColorOrder', lines)
        
        plot(bins, mean(sprate(goodTrials & stim.targchosen==1, :))); hold on
        plot(bins, mean(sprate(goodTrials & stim.targchosen==2, :)));
        xlim([-.1 1.3])
        xlabel('Time from Motion Onset')
        if kNeuron==1
            ylabel('Spike Rate')
        end
        % title('Choice Sorted PSTH', 'FontWeight', 'normal')
        % legend({'Choice 1', 'Choice 2'}, 'Location', 'BestOutside')
        
        % --- Get Whitened PTA
        
        % get residuals (subtract off the effect of motion onset)
        residuals = bsxfun(@minus, sprate, nanmean(sprate));
        Rs(goodTrials, :, kNeuron) = residuals(goodTrials,:);
        
%         % using raw rates, not residuals
%         Rs(goodTrials, :, kNeuron) = sprate(goodTrials,:);
        
        % pulse values
        pulses = sum(stim.pulses(goodTrials,:,:),3);
        % pulse times (in bins relative to spcnt column 1)
        pulseTimes = find(pdsa.binSpTimes(stim.timing(find(stim.goodtrial,1)).pulses - stim.timing(find(stim.goodtrial,1)).motionon, 0, window, binSize));
        
        nTimeLags = 60;
%        pta = pulseSTA(pulses, residuals(goodTrials,:), pulseTimes, nTimeLags);
        pta = pulseSTA(pulses, sprate(goodTrials,:), pulseTimes, nTimeLags);
        ptaBins = bsxfun(@plus, (1:nTimeLags)', pulseTimes) * binSize;
        
        
        subplot(3,nNeurons,2*nNeurons + kNeuron)
        
        set(gcf, 'DefaultAxesColorOrder', hot(12))
        plot(ptaBins, pta)
        axis tight
        % leg_str = regexp(sprintf('Pulse %d,', (1:7)'), ',', 'split');
        % legend(leg_str{1:7}, 'Location', 'BestOutside')
        xlabel('Time')
        if kNeuron ==1
            ylabel('Spikes / sec / gabor')
        end
        % title('PTA', 'FontWeight', 'normal')
    end % neuron loop
    drawnow
    set(gcf, 'PaperSize', [nNeurons*3 4], 'PaperPosition', [0 0 [nNeurons*3 4]])
    %saveas(gcf, [savePath filesep exname '.pdf'])
    
    nnPerSession(kEx) = nNeurons;
    
    %% Now do population decoding approach
    
    % Rs is trialxbinsxneuron - sum across time bins, and then squeeze.
    R = squeeze(sum(Rs(:,bins > 0 & bins < 1.2,:),2));
    Cho = sign(stim.targchosen - 1.5);
    Direc = sign(sum(sum(stim.pulses, 3),2));
    
    goodTrials = ~any(isnan(R),2);
    
    R = R(goodTrials,:);
    Cho = Cho(goodTrials);
    Direc = Direc(goodTrials);
    spikes = Rs(goodTrials,:,:);
    dirprob = stim.dirprob(goodTrials);
    
    if corr(Cho, Direc) < 0
        Cho = -Cho;
    end
    
    froIx = stim.trialId==stim.frozenTrialIds;
    froIx = froIx(goodTrials);    
    
    % get weights
    % use GLM fit
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
    end
    
    
    if corr(R*w2(2:end), Cho) < 0
        flipDV = true;
    else
        flipDV = false;
    end
    
    wAll = w2(2:end);
    kAll = w2(1);
        
    allW{kEx} = wAll;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % get a decision variable for each trial
%     [nTrials, nBins, nNeurons] = size(spikes);
%     
%     dvAll   = nan(nTrials, nBins);
% 
%     for iT = 1:nBins
%         R = squeeze(spikes(:,iT,:) );
%         dvAll(:,iT)   = R*wAll ;
% 
%         if flipDV
%             dvAll(:,iT)   = -dvAll(:,iT);
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get a decision variable for each trial
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
    
    %[cp_m,cp_s] = choiceProbabilityCalculate(dvAll(froIx, :), binaryCho(froIx));
    [cp_m,cp_s] = choiceProbabilityCalculate(dvAll, binaryCho(froIx));


       
    
    %% package and save
    
    sessionStruct.dvAll     = dvAll;
%     sessionStruct.dvTrain   = dvTrain;
%     sessionStruct.dvTest     = dvTest;

    sessionStruct.bins       = bins;
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
    sessionStruct.flipDV   = flipDV;
    sessionStruct.nNeurons = nNeurons;
    

%     if strcmp(condition{2}, 'early')
%         save([savePath(1:end-5) 'data' filesep 'allRevco' filesep exname], '-v7.3', '-struct', 'sessionStruct')
%     else
%         save([savePath(1:end-4) 'data' filesep 'allRevco' filesep exname], '-v7.3', '-struct', 'sessionStruct')
%     end
    

    if strcmp(condition{2}, 'early')
        save([savePath(1:end-5) 'data' filesep exname], '-v7.3', '-struct', 'sessionStruct')
    else
        save([savePath(1:end-4) 'data' filesep exname], '-v7.3', '-struct', 'sessionStruct')
    end
    
end % experiment loop




