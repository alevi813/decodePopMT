
function cpt_bySession(model, condition, fitWindow, plotWindow)

xvalPath = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/sessionData/' model filesep fitWindow 'Fit_' plotWindow 'Plot'];
allPath  = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/sessionData/no_xval/' model filesep fitWindow 'Fit_' plotWindow 'Plot'];
%%
experiments = getExperimentsAnd({condition});
nExp = length(experiments);

for iExp = 1:nExp
    
    xval_session = load([xvalPath filesep experiments{iExp} '.mat']);
    all_session = load([allPath filesep experiments{iExp} '.mat']);
    
    figure(1)
    subplot(3,ceil(nExp/3),iExp)
    plot(xval_session.bins, smooth(xval_session.cp_m, 10)); hold on
    plot(all_session.bins, smooth(all_session.cp_m, 10));
    title(experiments{iExp})
    if strcmp(plotWindow, 'stimulus')
        xlim([-.5 1.5])
    else
        xlim([-1 0])
    end
    if iExp ==1
        ylabel('CP')
        xlabel('time')
        
        legend('xval', 'no xval')
    end
    
    figure(2)
    subplot(3,ceil(nExp/3),iExp)
    imagesc(cell2mat(arrayfun(@(x) x.wTrain , xval_session.sLoo, 'UniformOutput', false))); hold on
    title(experiments{iExp})
    if iExp ==1
        xlabel('nTrial')
        ylabel('nNeuron')
    end
    
    figure(3)
    subplot(3,ceil(nExp/3),iExp)
    %histogram(cell2mat(arrayfun(@(x) x.wTrain , xval_session.sLoo, 'UniformOutput', false))); hold on
    plot(mean(cell2mat(arrayfun(@(x) x.wTrain , xval_session.sLoo, 'UniformOutput', false)), 2), '*'); hold on
    plot(all_session.wAll, 'o');
    title(experiments{iExp})
    if iExp ==1
        ylabel('weight')
        xlabel('n neuron')
        
        legend('xval (avg)', 'no xval')
    end    
end