% test validation

model       = 'direction';
condition   = 'flat';
experiments = getExperimentsAnd({condition});

dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model '/data'];
figPath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model filesep condition];

doPlot = true;

for iExp = 1:length(experiments)
    % load session
    session = load([dataPath filesep experiments{iExp} '.mat']);
    
    % get sumcoh
    if isfield(session.stim, 'validTrials')
        pulses = session.stim.pulses(session.stim.validTrials, :,:);
    else
        pulses = session.stim.pulses(session.stim.goodtrial, :,:);
    end
    pulses = sum(pulses, 3) ./19;
    sumcoh = mean(pulses,2);
    nBins  = 10;
    cohbins = binCoherences(sumcoh, nBins); % bin cohs for later
    
    % rename weight vectors
    trainData = session.wTrain;
    testData = session.wTest;
    
    % get pearson's r on two sets of weights
    if length(trainData) > 1
        tmpr = corrcoef(trainData, testData);
        rval(iExp) = tmpr(1,2);
    else
        rval(iExp) = NaN;
    end
    
    % Calculate predition accuracy for choice or direction using ROC. For
    % direction, calculate a full neurometric function and get the
    % threshold
    %if strcmp(model, 'direction')
        noDir = session.direction==0;
        direction = session.direction(~noDir) ==1;
%        wSpikes   = mean(session.dvTrain(~noDir, (31:102)), 2);
        wSpikes   = mean(session.dvAll(~noDir, (31:102)), 2);
        sumcoh    = sumcoh(~noDir);
        binId     = cohbins.id(~noDir);
        
        binsub = [0 1 2 3 4];
        
        for kBin = 1:nBins/2
            tmpdir    = [direction(binId==kBin); direction(binId==(nBins-binsub(kBin)))];
            tmpspikes = [wSpikes(binId==kBin); wSpikes(binId==(nBins-binsub(kBin)))];
            
            if length(unique(tmpdir)) > 1
%                [~, ~, ~, bin_auc(kBin)] = perfcurve(tmpdir, tmpspikes, 1);
                tmp(kBin) = roc(tmpspikes, tmpdir);
                bin_auc(kBin) = tmp(kBin).AUC;
                aucse(kBin) = tmp(kBin).serror;
            else
                bin_auc(kBin) = NaN;
            end
            
            if bin_auc(kBin) < .5
                bin_auc(kBin) = 1-bin_auc(kBin);
            end
        end
        
%         if iExp == 10
%             bin_auc(3) = NaN;
%         end
        
        % fit a logarithmic function to the auc points
        xvals = abs(cohbins.binCenters(1:nBins/2)); xfit = linspace(xvals(1),xvals(end));
        xvals(isnan(bin_auc)) = []; aucse(isnan(bin_auc)) = []; bin_auc(isnan(bin_auc)) = []; 
        logfun = @(x,xvals) (x(1)+x(2)*log(xvals)); x0 = [0, .5];
        xhat = lsqcurvefit(logfun,x0,xvals,bin_auc');
        fitvals = logfun(xhat,xfit);
        
        if doPlot
            if strcmp(condition, 'early')
                subplot(2, ceil(length(experiments)/2), iExp); hold on
            else
                subplot(1, length(experiments), iExp); hold on
            end
%            plot(xvals,bin_auc,'ko', xfit, fitvals,'b-');
            errorbar(xvals,bin_auc, aucse, 'k+'); 
            plot(xfit, fitvals,'r--');
            ylim([0.5 1.01])
            yticks([.5:.25:1])
            title(experiments{iExp})
            supertitle([condition '-stimulus session nmfs'], 12)
            %set(gcf, 'yticklabels', [0 .75 1])
            %axis('square')
        end
        
        
        % get threshold of your neurometric function
        [~, threshIx] = min(abs(fitvals-.75));
        nmThresh(iExp) = xfit(threshIx);
        %text(.25, .75, num2str(nmThresh(iExp)));
    %else
        %roc = roc(mean(session.dvTrain(session.froIx, (31:102)), 2), session.cho(session.froIx));
        [~, ~, ~, auc_tmp] = perfcurve(session.cho(session.froIx), mean(session.dvAll(session.froIx, (31:102)), 2), 1);
        
        if auc_tmp < .5
            cp(iExp) = 1-auc_tmp;
        else
            cp(iExp) = auc_tmp;
        end
    %end % model check
    
end % exp loop

if strcmp(model, 'direction')
    dir.(condition).wCorr = rval;
%    dir.(condition).auc   = auc;
    dir.(condition).nmThresh   = nmThresh;
    dir.(condition).cp   = cp;
else
    cho.(condition).wCorr = rval;
    cho.(condition).cp   = cp;
    cho.(condition).nmThresh   = nmThresh;    
end