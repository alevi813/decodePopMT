% test validation

comp   = getComp;
doPlot = false;

model       = 'choice';
condition   = {'flat', 'early', 'late'};
%subject     = 'leo';

for iCond = 1:length(condition)
    experiments = getExperimentsAnd(condition{iCond});
    
    %experiments = getExperimentsAnd( {subject, condition{iCond}} );
    
    if strcmp(comp, 'laptop')
        dataPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model '/rawRates/data'];
        figPath   = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/' model filesep condition{iCond}];
    else
        dataPath  = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/' model '/data'];
        figPath   = ['/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/' model filesep condition{iCond}];
    end
    
    for iExp = 1:length(experiments)
        nBins  = 10;
        % load session
        session = load([dataPath filesep experiments{iExp} '.mat']);
        
        if session.nNeurons > 1 %else do nothing
            % get sumcoh
            if isfield(session.stim, 'validTrials')
                pulses = session.stim.pulses(session.stim.validTrials, :,:);
            else
                pulses = session.stim.pulses(session.stim.goodtrial, :,:);
            end
            pulses = sum(pulses, 3) ./19;
            sumcoh = mean(pulses,2);
            sumcoh = zscore(sumcoh);
            
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
            
            % Calculate prediction accuracy for choice or direction using ROC. For
            % direction, calculate a full neurometric function and get the
            % threshold
            %if strcmp(model, 'direction')
            noDir = session.direction==0;
            direction = session.direction(~noDir) ==1;
            wSpikes   = mean(session.dvAll(~noDir, (50:170)), 2);        whichCond = condition{iCond};
            %         switch whichCond
            %             case 'flat'
            %                 wSpikes   = mean(session.dvAll(~noDir, (50:170)), 2);
            %             case 'late'
            %                 wSpikes   = mean(session.dvAll(~noDir, (110:170)), 2);
            %             case 'early'
            %                 wSpikes   = mean(session.dvAll(~noDir, (50:110)), 2);
            %         end
            
            sumcoh    = sumcoh(~noDir);
            binId     = cohbins.id(~noDir);
            
            %            binsub = [0 1 2 3 4 5];
            binsub = [0 1 2 3 4];
            
            for kBin = 1:nBins/2
                tmpdir    = [direction(binId==kBin); direction(binId==(nBins-binsub(kBin)))];
                tmpspikes = [wSpikes(binId==kBin); wSpikes(binId==(nBins-binsub(kBin)))];
                
                if length(unique(tmpdir)) > 1
                    %                [~, ~, ~, bin_auc(kBin)] = perfcurve(tmpdir, tmpspikes, 1);
                    tmp(kBin) = roc(tmpspikes, tmpdir);
                    bin_auc(kBin) = tmp(kBin).AUC;
                    %aucse(kBin) = tmp(kBin).serror;
                else
                    bin_auc(kBin) = NaN;
                end
                
                if bin_auc(kBin) < .5
                    bin_auc(kBin) = 1-bin_auc(kBin);
                end
            end
            
            % fold/unsign your cohbins
            fitvals = abs(cohbins.binCenters(1:nBins/2));
            % save single expt vals into larger matrices
            xAll(iExp, :) = fitvals;
            allAuc(iExp, :) = bin_auc;
            
            % get threshold of your neurometric function
            %     [~, threshIx] = min(abs(fitvals-.75));
            %     nmThresh(iExp) = xfit(threshIx);
            %text(.25, .75, num2str(nmThresh(iExp)));
            %else
            %roc = roc(mean(session.dvTrain(session.froIx, (31:102)), 2), session.cho(session.froIx));
            
            % dprime
            dp(iExp) = dprime(mean(session.newdvAll(:, (50:170)), 2), direction);
            %dp(iExp) = dprime(mean(session.dvAll(:, (50:170)), 2), direction);
            
            % whole trial cp
            [~, ~, ~, auc_tmp] = perfcurve(session.cho(session.froIx), mean(session.newdvAll(session.froIx, (50:170)), 2), 1);
            % first half
            %[~, ~, ~, auc_tmp] = perfcurve(session.cho(session.froIx), mean(session.newdvAll(session.froIx, (50:110)), 2), 1);
            % second half
            %[~, ~, ~, auc_tmp] = perfcurve(session.cho(session.froIx), mean(session.newdvAll(session.froIx, (110:170)), 2), 1);
            
            
            if auc_tmp < .5
                cp(iExp) = 1-auc_tmp;
            else
                cp(iExp) = auc_tmp;
            end
            %end % model check
            % all trial direction prediction accuracy.
            tmpS = roc(wSpikes, direction);
            accuracy(iExp) = tmpS.AUC;
        else
            xAll(iExp, :)   = nan(nBins/2,1);
            allAuc(iExp, :) = nan(nBins/2,1);
            cp(iExp)        = NaN;
            accuracy(iExp)  = NaN;
        end % if >1 neuron
        
    end % exp loop
    
    
    if strcmp(model, 'direction')
        %dir.(condition{iCond}).wCorr = rval;
        dir.(condition{iCond}).nmf.auc    = allAuc;
        dir.(condition{iCond}).nmf.xvals  = xAll;
        %dir.(condition{iCond}).nmThresh   = nmThresh;
        dir.(condition{iCond}).cp         = cp;
        dir.(condition{iCond}).dp         = dp;
        dir.(condition{iCond}).accuracy   = accuracy;
    else
        %cho.(condition{iCond}).wCorr = rval;
        cho.(condition{iCond}).cp   = cp;
        cho.(condition{iCond}).dp   = dp;
        cho.(condition{iCond}).nmf.auc = allAuc;
        cho.(condition{iCond}).nmf.xvals = xAll;
        %cho.(condition{1}).nmThresh   = nmThresh;
        cho.(condition{iCond}).accuracy   = accuracy;
    end
    
    clearvars -except dir cho comp doPlot model condition iCond subject nBins nancy leo
    
end %iCond
%%
if doPlot
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(232); hold on
    % flat /// direction
    dir.flat.nmf.m.x = nanmean(dir.flat.nmf.xvals);
    dir.flat.nmf.m.auc = nanmean(dir.flat.nmf.auc);
    dir.flat.nmf.se = nanstd(dir.flat.nmf.auc)/sqrt(size(dir.flat.nmf.auc, 1));
    
    % fit a logarithmic function to the auc points
    [xfit, fitvals] = logFit(dir.flat.nmf.m.x, dir.flat.nmf.m.auc);
    
    errorbar(dir.flat.nmf.m.x, dir.flat.nmf.m.auc, dir.flat.nmf.se, 'ko');
    plot(xfit, fitvals, 'color', [0 0 .8], 'linewidth', 1.5)
    
    % flat /// choice
    cho.flat.nmf.m.x = nanmean(cho.flat.nmf.xvals);
    cho.flat.nmf.m.auc = nanmean(cho.flat.nmf.auc);
    cho.flat.nmf.se = nanstd(cho.flat.nmf.auc)/sqrt(size(cho.flat.nmf.auc, 1));
    
    % fit a logarithmic function to the auc points
    [xfit, fitvals] = logFit(cho.flat.nmf.m.x, cho.flat.nmf.m.auc);
    
    errorbar(cho.flat.nmf.m.x, cho.flat.nmf.m.auc, cho.flat.nmf.se, 'k^');
    plot(xfit, fitvals, 'color', [0 0 .3], 'linewidth', 1.5)
    axis square
    xlabel('Motion strength (z)', 'fontsize', 14)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(233); hold on
    % late /// direction
    dir.late.nmf.m.x = nanmean(dir.late.nmf.xvals);
    dir.late.nmf.m.auc = nanmean(dir.late.nmf.auc);
    dir.late.nmf.se = nanstd(dir.late.nmf.auc)/sqrt(size(dir.late.nmf.auc, 1));
    
    % fit a logarithmic function to the auc points
    [xfit, fitvals] = logFit(dir.late.nmf.m.x, dir.late.nmf.m.auc);
    
    errorbar(dir.late.nmf.m.x, dir.late.nmf.m.auc, dir.late.nmf.se, 'ko');
    plot(xfit, fitvals, 'color', [.8 .8 0], 'LineWidth', 1.5)
    
    % late /// choice
    cho.late.nmf.m.x = nanmean(cho.late.nmf.xvals);
    cho.late.nmf.m.auc = nanmean(cho.late.nmf.auc);
    cho.late.nmf.se = nanstd(cho.late.nmf.xvals)/sqrt(size(cho.late.nmf.xvals, 1));
    
    % fit a logarithmic function to the auc points
    [xfit, fitvals] = logFit(cho.late.nmf.m.x, cho.late.nmf.m.auc);
    
    errorbar(cho.late.nmf.m.x, cho.late.nmf.m.auc, cho.late.nmf.se, 'k^');
    plot(xfit, fitvals, 'color', [.4 .4 0], 'LineWidth', 1.5)
    axis square
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(231); hold on
    % early /// direction
    dir.early.nmf.m.x = nanmean(dir.early.nmf.xvals);
    dir.early.nmf.m.auc = nanmean(dir.early.nmf.auc);
    dir.early.nmf.se = nanstd(dir.early.nmf.auc)/sqrt(size(dir.early.nmf.auc, 1));
    
    % fit a logarithmic function to the auc points
    [xfit, fitvals] = logFit(dir.early.nmf.m.x, dir.early.nmf.m.auc);
    
    errorbar(dir.early.nmf.m.x, dir.early.nmf.m.auc, dir.early.nmf.se, 'ko');
    plot(xfit, fitvals, 'color', [.8 0 0], 'linewidth', 1.5)
    
    % early /// choice
    cho.early.nmf.m.x = nanmean(cho.early.nmf.xvals);
    cho.early.nmf.m.auc = nanmean(cho.early.nmf.auc);
    cho.early.nmf.se = nanstd(cho.early.nmf.auc)/sqrt(size(cho.early.nmf.auc, 1));
    
    % fit a logarithmic function to the auc points
    [xfit, fitvals] = logFit(cho.early.nmf.m.x, cho.early.nmf.m.auc);
    
    errorbar(cho.early.nmf.m.x, cho.early.nmf.m.auc, cho.early.nmf.se, 'k^');
    plot(xfit, fitvals, 'color', [.3 0 0], 'linewidth', 1.5)
    axis square
    ylabel('Accuracy', 'fontsize', 14)
    
end