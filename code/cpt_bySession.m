
function cpt_bySession(model, condition, fitWindow, plotWindow)

dataPath = ['/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/sessionData/' model filesep fitWindow 'Fit_' plotWindow 'Plot'];
%%
experiments = getExperimentsAnd({condition});
nExp = length(experiments);

for iExp = 1:nExp
    
    session = load([dataPath filesep experiments{iExp} '.mat']);
    
%     figure(1)
%     subplot(3,ceil(nExp/3),iExp)
%     plot(session.bins, smooth(session.cp_m, 10)); hold on
%     title(experiments{iExp})
%     if strcmp(plotWindow, 'stimulus')
%         xlim([-.5 1.5])
%     else
%         xlim([-1 0])
%     end
%     if iExp ==1
%         ylabel('CP')
%         xlabel('time')
%     end
%     
%     figure(2)
%     subplot(3,ceil(nExp/3),iExp)
%     imagesc(cell2mat(arrayfun(@(x) x.wTrain , session.sLoo, 'UniformOutput', false))); hold on
%     title(experiments{iExp})
%     if iExp ==1
%         xlabel('nTrial')
%         ylabel('nNeuron')
%     end
%     
%     figure(3)
    subplot(3,ceil(nExp/3),iExp)
    plot(mean(cell2mat(arrayfun(@(x) x.wTrain , session.sLoo, 'UniformOutput', false)), 2), 'o', 'markersize', 8, 'linewidth',2); hold on
    plot(median(cell2mat(arrayfun(@(x) x.wTrain , session.sLoo, 'UniformOutput', false)), 2), '*'); 
    title(experiments{iExp})
    if iExp ==1
        ylabel('weight')
        xlabel('n neuron')
        
        %legend('mean', 'median')
    end    
end