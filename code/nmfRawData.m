% get avg PSTH for each condition,
%dataPath = getpref('mtlipglm', 'dataPath');
dataPath =     '/Users/aaronlevi/Dropbox/twagAnalysis4.1/Data/nancy';
condition = {'early'}; % flat/late/early
fitDir    = 'fits_coupling2'; % main_fits / fits_notargs / fits_pulsecovars / fits_cho1 / fits_coupling / fits_noHist

S = loadUpFits(dataPath, 'MT', condition, fitDir);

%% Plot Model PSTH on top of data PSTH
% The three models compared in figure 5 are the "Stimulus-to-LIP" model,
% the "MT-to-LIP" model and "MT-to-LIP (full choice)"

modelNames=arrayfun(@(x) x.name, S(1).model, 'UniformOutput', false);
dp = arrayfun(@(x) x.model(1).dprime, S);
dpix = abs(dp) >= 0.1; %0.15

% The three models have different names in the code:
% Stimulus-to-LIP = Poisson
% MT-to-LIP       = MTsim
% MT-to-LIP (full choice) = MTsimChoice
% Compare them to the data in each plot
modelIxs={[find(strcmp(modelNames, 'data'))]};

S = S(dpix);
unitExpName = arrayfun(@(x) x.exname, S, 'UniformOutput', false);
expName = unique(unitExpName);
dp = dp(dpix);

for iModel=1:numel(modelIxs)
    %figure(iModel); clf
    
    modelIx  = modelIxs{iModel};
    nModels  = numel(modelIx);
    nNeurons = numel(S);
    
    plotFields={'psthCoh'};
    nPlotFields = numel(plotFields);
    
    for f=1:nPlotFields
        %        subplot(1,nPlotFields, f)
        
        field = 'trialSpikes';
        % --- Change Color Map
        if any(strfind(field, 'pta'))
            cmap=hot(12);
        else
            cmap=cbrewer('jake', 'rdbu', 8);
        end
        
        for kM=1:nModels
            
            kModel=modelIx(kM);
            
            % get trial spike trains for each neuron
            trialSpikes = arrayfun(@(x) x.model(kModel).(field), S, 'UniformOutput', false);
            timex=S(1).model(1).psthTime(:); stimOn = 101:211;
            
            for kN = 1:nNeurons
               nTrials(kN) = length(trialSpikes{kN});
            end
            
            for kSession = 1:length(expName)
                % get stim file, pulses
                stim    = getStim(expName{kSession}, dataPath);
                pulses  = mean(stim.pulses, 3);
                pulses  = mean(pulses, 2);
                if isfield(stim, 'validTrials')
                    pulses = pulses(stim.validTrials);
                else
                    pulses  = pulses(stim.goodtrial);
                end
                
                isFro   = pulses == 0;                
                
                % get units organized by session
                expIx = find(strcmp(unitExpName, expName{kSession}));                
                for iExp = 1:length(expIx)
                    % tmpSpikes{iExp} = trialSpikes{expIx(iExp)};
                    % tmpSpikes = cell2mat(tmpSpikes');
                    tmpSpikes = trialSpikes{expIx(iExp)};

                     % give signed direction according to pref dir (from dp)
                    if dp(iExp) < 1
                        direction = pulses < 0;
                    else 
                        direction = pulses > 0;
                    end
                    direction(isFro) = [];                
                    tmpSpikes(isFro, :) = [];
                    
                    auc{kSession}(iExp) = roc(mean(tmpSpikes(:,stimOn), 2), direction');
                    
                    if auc{kSession}(iExp).AUC < .5
                        auc{kSession}(iExp).AUC = 1 - auc{kSession}(iExp).AUC;
                    end
                end % iExp
                
                
            end % kSession
%             % normalize by max
%             cohpData = reshape(cell2mat(arrayfun(@(x) x.model(kModel).(field)/max(x.model(1).(field)(:)), S, 'UniformOutput', false)), [size(S(1).model(kModel).(field)) nNeurons]);
%             
%             
%             % cohpData = cohpData(:,:,dpix);
%             %popcoh = nanmean(cohpData, 3);
%             popcohDir1   = nanmean( nanmean(cohpData(:,1:4, :),3), 2 );
%             popcohDir2   = nanmean( nanmean(cohpData(:,5:8, :),3), 2 );
%             
%             
%             timex=S(1).model(1).psthTime(:); stimOn = 101:211;
%             xdim = [-500 1500];
%             
%             binsub = [0 1 2 3 4];
%             nBins = length(S(1).model(1).CohEdges);
%             for kBin = 1:nBins/2
%                 tmpspikes = [popcoh(stimOn, kBin), popcoh(stimOn, (nBins-binsub(kBin)))];
%                 %tmpdir    = [ones(length(direction(binId==(nBins-binsub(kBin)))];
%                 
%                 if length(unique(tmpdir)) > 1
%                     [~, ~, ~, bin_auc(kBin)] = perfcurve(tmpdir, tmpspikes, 1);
%                 else
%                     bin_auc(kBin) = NaN;
%                 end
%                 
%                 if bin_auc(kBin) < .5
%                     bin_auc(kBin) = 1-bin_auc(kBin);
%                 end
%             end
            
        end %nmodels
    end %nplotfields
    
end % imodel (modelixs)

%%

for iSession = 1:length(expName)   
    rawAccuracy.m(iSession)  = mean(cell2mat(arrayfun(@(x) x.AUC, auc{iSession}, 'UniformOutput', false))) ;
    rawAccuracy.se(iSession) = std(cell2mat(arrayfun(@(x) x.AUC, auc{iSession}, 'UniformOutput', false))) / sqrt(numel(  cell2mat(arrayfun(@(x) x.AUC, auc{iSession}, 'UniformOutput', false))) ); 
end
