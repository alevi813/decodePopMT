function dCond = decodingSimilarity(dCond)
% measure similarity between decoding models

%% similarity over time

% get cosine similarity
% dot product of the two model dvs, normalized by euclidean norm
thetaTime   = nan(1,124);
%     thetaTimeTrain = nan(1,120);
%     thetaTimeTest  = nan(1,120);

for iBin = 1:length(dCond.dir.bins)
    
    % all data ACROSS dir/cho
    tempDot     = dot(dCond.dir.newdv(:, iBin), dCond.cho.newdv(:, iBin));
    tempNormDir = norm(dCond.dir.newdv(:,iBin));
    tempNormCho = norm(dCond.cho.newdv(:,iBin));
    
    thetaTime(iBin) = acosd( tempDot / (tempNormDir*tempNormCho) );
    
    %         % WITHIN dir model
    %         tempDot       = dot(dvDirTrain(iBin,:), dvDirTest(iBin,:));
    %         tempNormTrain = norm(dvDirTrain(iBin,:));
    %         tempNormTest  = norm(dvDirTest(iBin,:));
    %
    %         thetaTime.within(iBin) = acosd( tempDot / (tempNormDir*tempNormTest) );
    
end % bin loop

dCond.thetaTime = thetaTime;

%%
% 
% doPlot = true;
% 
% if doPlot
%     figure
%     subplot(2,3,1); hold on
%     plot(session.bins, mean(combo.dvAll(combo.cho==1,:)))
%     plot(session.bins, mean(combo.dvAll(combo.cho==0,:)) .*-1 )
%     ylabel('weighted spks')
%     xlabel('Time')
%     xlim([-.5 1.6])
%     title('choice sorted psth')
%     
%     subplot(2,3,2); hold on
%     set(gcf, 'DefaultAxesColorOrder', hot(12))
%     plot(session.ptaBins, combo.mPTA, 'LineWidth', 1.75, 'Color', [.2 .2 .2])
%     title('PTA')
%     
%     subplot(2,3,4); hold on
%     plot(session.bins, smooth(mean(combo.cpt), 9))
%     plot([-.5 1.5], [.5 .5], 'r--')
%     title('CP')
%     ylabel('CP')
%     xlabel('Time')
%     xlim([-.5 1.6])
%     title('choice probability')
%     
%     subplot(2,3,5)
%     plot(sum(combo.mPTA), 'k-o')
%     %ylabel('PTA pulse sum')
%     set(gca, 'Xtick', [1:7])
%     xlim([0 8])
%     title('PTA pulse sum')
%     
%     subplot(2,3,6)
%     errorbar(1:7, mean(combo.ppk), std(combo.ppk)/sqrt(length(experiments)), '-o' )
%     xlim([0 8])
%     set(gca, 'Xtick', [1:7])
%     title('PPK')
%     ylabel('weight')
%     xlabel('pulse')
%     
%     supertitle([model ' - ' condition ' - avg '], 12)
% end
% % saveas(gcf, [figPath filesep 'avg_popDV.pdf']);
% %%
% 
% % get some min/max info for the fancy rectangle
% mx(1) = max(cpDir.m);
% mn(1) = min(cpDir.m);
% mx(2) = max(cpChoice.m);
% mn(2) = min(cpChoice.m);
% 
% mx = max(mx); mn = min(mn);
% toplim    = mx + 0.005;
% botlim    = mn - 0.005;
% 
% % plot
% %figure; hold on
% rectangle('Position', [0, mn-0.005, 1.05, toplim- botlim], 'EdgeColor', [0.9 0.9 0.9], 'FaceColor', [0.9 0.9 0.9]);
% plot([session.bins(1) session.bins(end)], [0.5 0.5], 'k--')
% 
% boundedline(session.bins, cpDir.m, cpDir.se, 'cmap', [.8 0 0], 'alpha');
% boundedline(session.bins, cpChoice.m, cpChoice.se, 'cmap', [.4 0 0], 'alpha');
% 
% xlim([-.5 1.5])
% 
% %
% % late
% % [.8 .8 0]
% % [.4 .4 0]
