
function [sessionStruct] = sessionDecodingBackwards(model, condition, weightOrigin)
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
if nargin <3
    weightOrigin = 'calculateW';
end

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
    
    switch weightOrigin
        case 'calculateW'
            if strcmp(model, 'choice')
                savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/xval_rawRates/' condition{2}];
                %savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/backwardsWindow/' condition{2}];
                %savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/bw_fullWindow/' condition{2}];
            else
                savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/backwardsWindow/' condition{2}];
                %savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/bw_fullWindow/' condition{2}];
            end
        case 'loadW'
            if strcmp(model, 'choice')
                savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/bw_oldWeights/' condition{2}];
            else
                savePath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/bw_oldWeights/' condition{2}];
            end
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
    %window  = [-1 0]; % backwards windown when aligned to 'go'
    
    nBins = diff(window)/binSize;
    
    % spike rate for the population
    %Rs = nan(stim.nTrials, nBins, nNeurons); %preallocate
    Rs = nan(sum(stim.goodtrial), nBins, nNeurons); %preallocate
    
    % stimulus timing variables
    motionOnset  = [stim.timing(:).motionon] + [stim.timing(:).plxstart];
    motionOffset = [stim.timing(:).motionoff] + [stim.timing(:).plxstart];
    goTime       = [stim.timing(:).fpoff] + [stim.timing(:).plxstart];
    
    timeToGo = goTime - motionOffset;
    
    % loop over neurons in the session
    for kNeuron = 1:nNeurons
        
        
        spikeTimes   = neurons(kNeuron).spikeTimes;
        
        % get spike count aligned to motion onset
        %[spcnt, bins]  = pdsa.binSpTimes(spikeTimes, motionOnset, window, binSize);
        [spcnt, bins]  = pdsa.binSpTimes(spikeTimes, motionOnset, window, binSize);
        
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
        %xlim([-.1 1.3])
        xlabel('Time from Motion Onset')
        if kNeuron==1
            ylabel('Spike Rate')
        end
        % title('Choice Sorted PSTH', 'FontWeight', 'normal')
        % legend({'Choice 1', 'Choice 2'}, 'Location', 'BestOutside')
        
        
        % get residuals (subtract off the effect of motion onset)
        residuals = bsxfun(@minus, sprate, nanmean(sprate));
        Rs(goodTrials, :, kNeuron) = residuals(goodTrials,:);
        %Rs(:, :, kNeuron) = residuals(goodTrials,:);
        
        %         % using raw rates, not residuals
        %         Rs(goodTrials, :, kNeuron) = sprate(goodTrials,:);
        
        % pulse values
        pulses = sum(stim.pulses(goodTrials,:,:),3);
        
        
        subplot(3,nNeurons,2*nNeurons + kNeuron)
        
        %%%% REDO TEMPORAL WINDOW ON SPIKE RATES FOR DV
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [spcnt_fromGo, bins_fromGo]  = pdsa.binSpTimes(spikeTimes, goTime, [-1 0], binSize);
        
        % smooth spike count with a boxcar filter
        sm = 5; % size of boxcar (number of bins)
        tmp_fromGo = filter(boxcar(sm)/sm, 1, spcnt_fromGo');
        tmp_fromGo = flipud(tmp_fromGo);
        tmp_fromGo = filter(boxcar(sm)/sm, 1, tmp_fromGo);
        tmp_fromGo = flipud(tmp_fromGo);
        
        spcnt_fromGo = tmp_fromGo';
        
        sprate_fromGo = spcnt_fromGo/binSize;
        
%         % residuals
%         residuals_fromGo                  = bsxfun(@minus, sprate_fromGo, nanmean(sprate_fromGo));
%         Rs_fromGo(goodTrials, :, kNeuron) = residuals_fromGo(goodTrials,:);
%         spikes_fromGo                     = Rs_fromGo(goodTrials,:,:);

        % using raw rates, not residuals
        Rs_fromGo(goodTrials, :, kNeuron) = sprate_fromGo(goodTrials,:);
        spikes_fromGo                     = Rs_fromGo(goodTrials,:,:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end % neuron loop
    drawnow
    set(gcf, 'PaperSize', [nNeurons*3 4], 'PaperPosition', [0 0 [nNeurons*3 4]])
    %saveas(gcf, [savePath filesep exname '.pdf'])
    
    nnPerSession(kEx) = nNeurons;
    
    %% Now do population decoding approach
    
    % Rs is trialxbinsxneuron - sum across time bins, and then squeeze.
    R = squeeze(sum(Rs_fromGo(:,bins_fromGo > -.5 & bins_fromGo < 0,:),2)); % summed spikes over 500ms pre "go"
    %R = squeeze(sum(Rs_fromGo,2)); % summed spikes based on entire window (1 second before "go" -- double check binsptimes window though)
    Cho = sign(stim.targchosen - 1.5);
    Direc = sign(sum(sum(stim.pulses, 3),2));
    
    %goodTrials = ~any(isnan(R),2);
    
    R = R(goodTrials,:);
    Cho = Cho(goodTrials);
    Direc = Direc(goodTrials);
    
    dirprob = stim.dirprob(goodTrials);
    
    if corr(Cho, Direc) < 0
        Cho = -Cho;
    end
    
    froIx = stim.trialId==stim.frozenTrialIds;
    froIx = froIx(goodTrials);
    
    % get weights
    % use GLM fit
    switch weightOrigin
        case 'calculateW'
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
        case 'loadW'
            if strcmp(model, 'choice')
                tmpStruct = load(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/residuals/data/' exname '.mat']);
            else
                tmpStruct = load(['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/residuals/data/' exname '.mat']);
            end % if choice or dir
            
            wAll = tmpStruct.wAll;
            sLoo = tmpStruct.sLoo;
    end % switch weightOrigin
        
    if corr(R*wAll, Cho) < 0
        flipDV = true;
    else
        flipDV = false;
    end    
    
    %allW{kEx} = wAll;
    
    % get a decision variable for each trial
    [nTrials, nBins_fromGo, nNeurons] = size(spikes_fromGo);
    
    if strcmp(model, 'choice')
        dvAll   = nan(sum(froIx), nBins_fromGo);
    else
        dvAll   = nan(nTrials, nBins_fromGo);
   end
     
   % reweighted decision variable
   if strcmp(model, 'choice')
       froSpikesFromGo = spikes_fromGo(froIx,:,:);
       for iTrial = 1:sum(froIx)
           R = squeeze(froSpikesFromGo(iTrial,:,:) ); %squeeze to get spikes per bin per neuron during a specific trial
           %dvAll(iTrial,:)   = R*wAll; % muliply spike matrix (nBin x nNeuron) by weight vector (nNeurons x 1)
           dvAll(iTrial,:)   = R* sLoo(iTrial).wTrain; % muliply spike matrix (nBin x nNeuron) by weight vector (nNeurons x 1)
           
           if flipDV
               dvAll(iTrial,:)   = -dvAll(iTrial,:);
           end
       end% ntrials        
   else
       for iTrial = 1:nTrials
           R = squeeze(spikes_fromGo(iTrial,:,:) ); %squeeze to get spikes per bin per neuron during a specific trial
           dvAll(iTrial,:)   = R*wAll; % muliply spike matrix (nBin x nNeuron) by weight vector (nNeurons x 1)
           %dvAll(iTrial,:)   = R* sLoo(iTrial).wTrain; % muliply spike matrix (nBin x nNeuron) by weight vector (nNeurons x 1)
           
           if flipDV
               dvAll(iTrial,:)   = -dvAll(iTrial,:);
           end
       end% ntrials       
   end % choice v direction model
    
    % smooth the dv
    dvAll   = filter(ones(5,1)/5, 1, dvAll')';
    %     dvTrain = filter(ones(5,1)/5, 1, dvTrain')';
    %     dvTest  = filter(ones(5,1)/5, 1, dvTest')';
    
    % calculate traditional cp
    % Cho is already aligned to pref direction of neuron above, but you
    % still need to binarize
    binaryCho = Cho==1;
    
    if strcmp(model, 'choice')
        [cp_m,cp_s] = choiceProbabilityCalculate(dvAll, binaryCho(froIx));
    else
        [cp_m,cp_s] = choiceProbabilityCalculate(dvAll(froIx, :), binaryCho(froIx));
    end
       
    plotOn = false;
    if plotOn
        figure(2); clf
        set(gcf, 'DefaultAxesColorOrder', lines)
        subplot(2,3,1)
        plot(bins_fromGo, mean(dvAll(Cho==1,:))); hold on
        plot(bins_fromGo, mean(dvAll(Cho==-1,:)) .* -1 );
        %    plot(bins, mean(dvAll(Cho==-1,:)));
        
        ylabel('DV')
        xlabel('Time')
        axis tight
        %xlim([-.5 1.6])
        title('Choice Sorted DV')
        
        subplot(2,3,2)
        
        % motionOn = repmat([1 zeros(1, numel(bins)-1)], size(dv,1),1);
        % motionOn = find(motionOn(:));
        % coh = sum(pulses,2);
        % cho = Cho==1;
        
        
        subplot(2,3,6)
        plot(S.pk.mle(1:7), 'k-o'); hold on
        plot(S.pk.mle(1:7)+S.pk.sd(1:7), '--k');
        plot(S.pk.mle(1:7)-S.pk.sd(1:7), '--k');
        axis tight
        title('PPK')
        set(gca, 'Xtick', [1:7])
        
        % % try plotting the population PTA effect
        % plot(sum(pta)*20, 'b-o')
        
        % subplot(2,3,4)
        % plot(pulse_dvR, 'b-o'); hold on
        % plot(pulse_dvL, 'r-o'); hold on
        
        subplot(2,3,4)
        plot(bins_fromGo, cp_m);
        %xlim([-.5 1.5])
        hold on
        %plot([-.5 1.5], [.5 .5], 'r--')
        title('CP')
        ylabel('CP')
        xlabel('Time')
        
        
        if strcmp(model, 'choice')
            supertitle(['Choice - ' condition], 12)
        else
            supertitle(['Direction - ' condition], 12)
        end
        %saveas(gcf, [savePath filesep exname '_popDV.pdf'])
        
    end %plotOn?
    %% package and save
    
    sessionStruct.dvAll      = dvAll;
    sessionStruct.bins       = bins;
    sessionStruct.bins_fromGo       = bins_fromGo;
    sessionStruct.spikes_fromGo     = spikes_fromGo;
    %sessionStruct.newdvAll          = [dvAll(binaryCho==1,:); (dvAll(binaryCho==0,:) .*-1)];
    %    sessionStruct.normdv    = mean(sessionStruct.newdv) / max(mean(sessionStruct.newdv));
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
    sessionStruct.timeToGo  = timeToGo;
    
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
    
    clearvars -except model condition weightOrigin savePath experiments comp nFiles nFilesNames nDataPath dataPath
end % experiment loop




