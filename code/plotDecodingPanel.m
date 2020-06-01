%% psth
figure(3)
subplot(232)
% plot(dS.flat.dir.bins, mean(dS.flat.cho.newdv), 'color', [0 0 .4]); hold on
% plot(dS.flat.dir.bins, mean(dS.flat.dir.newdv), 'color', [0 0 .8]);

dS.flat.cho.normdv = mean(dS.flat.cho.newdv) ./ max(mean(dS.flat.dir.newdv));
dS.flat.dir.normdv = mean(dS.flat.dir.newdv) ./ max(mean(dS.flat.dir.newdv));

plot(dS.flat.dir.bins, dS.flat.cho.normdv, 'color', [0 0 .4]); hold on
plot(dS.flat.dir.bins, dS.flat.dir.normdv, 'color', [0 0 .8]);
xlim([-.5 1.5])
axis square
%legend('choice', 'direction')

subplot(233)
% plot(dS.late.dir.bins, mean(dS.late.cho.newdv), 'color', [.4 .4 0]); hold on
% plot(dS.late.dir.bins, mean(dS.late.dir.newdv), 'color', [.8 .8 0]);

dS.late.cho.normdv = mean(dS.late.cho.newdv) ./ max(mean(dS.late.dir.newdv));
dS.late.dir.normdv = mean(dS.late.dir.newdv) ./ max(mean(dS.late.dir.newdv));

plot(dS.late.dir.bins, dS.late.cho.normdv, 'color', [.4 .4 0]); hold on
plot(dS.late.dir.bins, dS.late.dir.normdv, 'color', [.8 .8 0]);
xlim([-.5 1.5])
axis square
%legend('choice', 'direction')

subplot(231)
% plot(dS.early.dir.bins, mean(dS.early.cho.newdv), 'color', [.4 0 0]); hold on
% plot(dS.early.dir.bins, mean(dS.early.dir.newdv), 'color', [.8 0 0]);

dS.early.cho.normdv = mean(dS.early.cho.newdv) ./ max(mean(dS.early.dir.newdv));
dS.early.dir.normdv = mean(dS.early.dir.newdv) ./ max(mean(dS.early.dir.newdv));

plot(dS.early.dir.bins, dS.early.cho.normdv, 'color', [.4 0 0]); hold on
plot(dS.early.dir.bins, dS.early.dir.normdv, 'color', [.8 0 0]);
xlim([-.5 1.5])
axis square
%ylim([-1e-3 3.5e-3])
%legend('choice', 'direction')

set(gcf, 'color', 'w')

%% cpt
%figure

figure
%subplot(234); hold on
subplot(231); hold on
dS.early.dir.cpm  = smooth( mean(dS.early.dir.cpt), 10);
dS.early.dir.cpse = smooth( (std(dS.early.dir.cpt) / sqrt(size(dS.early.cho.cpt, 1))), 9);
boundedline(session.bins, dS.early.dir.cpm, dS.early.dir.cpse, 'cmap', [.8 0 0], 'alpha');
%boundedline(session.bins_fromGo, dS.early.dir.cpm, dS.early.dir.cpse, 'cmap', [.8 0 0], 'alpha');

dS.early.cho.cpm  = smooth( mean(dS.early.cho.cpt), 10);
dS.early.cho.cpse = smooth( (std(dS.early.cho.cpt) / sqrt(size(dS.early.cho.cpt, 1))), 9);
boundedline(session.bins, dS.early.cho.cpm, dS.early.cho.cpse, 'cmap', [.3 0 0], 'alpha');
%boundedline(session.bins_fromGo, dS.early.cho.cpm, dS.early.cho.cpse, 'cmap', [.3 0 0], 'alpha');

%xlim([-.5 1.4])
%ylim([0.46 .65])
axis square

%subplot(235); hold on
subplot(232); hold on
dS.flat.dir.cpm  = smooth( mean(dS.flat.dir.cpt), 10);
dS.flat.dir.cpse = smooth( (std(dS.flat.dir.cpt) / sqrt(size(dS.flat.cho.cpt, 1))), 9);
boundedline(session.bins, dS.flat.dir.cpm, dS.flat.dir.cpse, 'cmap', [0 0 .8], 'alpha');
%boundedline(session.bins_fromGo, dS.flat.dir.cpm, dS.flat.dir.cpse, 'cmap', [0 0 .8], 'alpha');

dS.flat.cho.cpm  = smooth( mean(dS.flat.cho.cpt), 10);
dS.flat.cho.cpse = smooth( (std(dS.flat.cho.cpt) / sqrt(size(dS.flat.cho.cpt, 1))), 9);
boundedline(session.bins, dS.flat.cho.cpm, dS.flat.cho.cpse, 'cmap', [0 0 .3], 'alpha');
%boundedline(session.bins_fromGo, dS.flat.cho.cpm, dS.flat.cho.cpse, 'cmap', [0 0 .3], 'alpha');

%xlim([-.5 1.5])
%ylim([0.46 .65])
axis square

%subplot(236); hold on
subplot(233); hold on
dS.late.dir.cpm  = smooth( mean(dS.late.dir.cpt), 10);
dS.late.dir.cpse = smooth( (std(dS.late.dir.cpt) / sqrt(size(dS.late.cho.cpt, 1))), 9);
boundedline(session.bins, dS.late.dir.cpm, dS.late.dir.cpse, 'cmap', [.8 .8 0], 'alpha');
%boundedline(session.bins_fromGo, dS.late.dir.cpm, dS.late.dir.cpse, 'cmap', [.8 .8 0], 'alpha');

dS.late.cho.cpm  = smooth( mean(dS.late.cho.cpt), 10);
dS.late.cho.cpse = smooth( (std(dS.late.cho.cpt) / sqrt(size(dS.late.cho.cpt, 1))), 9);
boundedline(session.bins, dS.late.cho.cpm, dS.late.cho.cpse, 'cmap', [.4 .4 0], 'alpha');
%boundedline(session.bins_fromGo, dS.late.cho.cpm, dS.late.cho.cpse, 'cmap', [.4 .4 0], 'alpha');

%xlim([-.5 1.5])
%ylim([0.46 .65])
axis square

%%%% for stats
dirall    = [dS.flat.dir.cpt; dS.late.dir.cpt; dS.early.dir.cpt];
dirallm   = smooth( mean(dirall), 10);
dirGrandM = mean(dirallm(50:170));
dirse = smooth( (std([dS.late.dir.cpt; dS.flat.dir.cpt; dS.early.dir.cpt]) / sqrt(size([dS.late.dir.cpt; dS.flat.dir.cpt; dS.early.dir.cpt], 1))), 9);
dirci = mean(dirse)*1.96;

choall    = [dS.flat.cho.cpt; dS.late.cho.cpt; dS.early.cho.cpt];
choallm   = smooth( mean(choall), 10);
choGrandM = mean(choallm(50:170));
chose = smooth( (std([dS.late.cho.cpt; dS.flat.cho.cpt; dS.early.cho.cpt]) / sqrt(size([dS.late.cho.cpt; dS.flat.cho.cpt; dS.early.cho.cpt], 1))), 9);
choci = mean(chose)*1.96;
%% alignment / similarity

% figure(1)
% subplot(234)
% plot(early.dir.bins, early.thetaTime, 'color', [.8 0 0])
% xlim([-.5 1.5])
% axis square
% 
% subplot(235)
% plot(flat.dir.bins, flat.thetaTime, 'color', [0 0 .8])
% xlim([-.5 1.5])
% axis square
% 
% subplot(236)
% plot(late.dir.bins, late.thetaTime, 'color', [.6 .8 0])
% xlim([-.5 1.5])
% axis square

%% pta

figure(3)
subplot(234); hold on
%subplot(131); hold on
plot(session.ptaBins, dS.early.dir.mPTA, 'LineWidth', 1.75, 'Color', [.8 0 0])
plot(session.ptaBins, dS.early.cho.mPTA, 'LineWidth', 1.75, 'Color', [.4 0 0])
axis square

subplot(235); hold on
%subplot(132); hold on
plot(session.ptaBins, dS.flat.dir.mPTA, 'LineWidth', 1.75, 'Color', [0 0 .9])
plot(session.ptaBins, dS.flat.cho.mPTA, 'LineWidth', 1.75, 'Color', [0 0 .4])
axis square

subplot(236); hold on
%subplot(133); hold on
plot(session.ptaBins, dS.late.dir.mPTA, 'LineWidth', 1.75, 'Color', [.8 .8 0])
plot(session.ptaBins, dS.late.cho.mPTA, 'LineWidth', 1.75, 'Color', [.4 .4 0])
axis square

set(gcf, 'color', 'white');

%% raw data scatter plots
% 
% figure
% % pulse cp vs. ppk
% 
% % pta peaks vs. ppk
% 
% %% direction model scatter plots
% 
% figure(1)
% figure(2)
% % pulse cp vs. ppk
% % get avg cp during a 150ms pulse
% for cond = 1:3
%     
%     figure(1)
%     if cond == 1
%         tmpcp  = flat.dir.cpm;        
%         ppk.m  = mean(flat.dir.combo.ppk);
%         ppk.se = std(flat.dir.combo.ppk) / sqrt(size(flat.dir.combo.ppk, 1));
%         
%         subplot(1,3,2); hold on
%         clr = winter(7);
%         title('Flat', 'FontSize', 15)
%     end
%     
%     if cond == 2
%         tmpcp = late.dir.cpm;
%         ppk.m = mean(late.dir.combo.ppk);
%         ppk.se = std(late.dir.combo.ppk) / sqrt(size(flat.dir.combo.ppk, 1));
%         
%         subplot(1,3,3); hold on
%         clr = autumn(7);
%         title('Late', 'FontSize', 15)
%     end
%     
%     if cond == 3
%         tmpcp = early.dir.cpm;
%         ppk.m = mean(early.dir.combo.ppk);
%         ppk.se = std(early.dir.combo.ppk) / sqrt(size(flat.dir.combo.ppk, 1));
%         
%         subplot(1,3,1); hold on
%         clr = hot(12);
%         title('Early', 'FontSize', 15)
%     end
%     
%     addPulseTime = [0 .150 .300 .450 .600 .750 .900];
%     for iPulse = 1:7
%         
%         cpPulse.m(iPulse) = mean(tmpcp(session.bins >=(0.05+addPulseTime(iPulse)) & session.bins <= (0.2+addPulseTime(iPulse))));
%         cpPulse.se(iPulse) = std(tmpcp(session.bins >=(0.05+addPulseTime(iPulse)) & session.bins <= (0.2+addPulseTime(iPulse)))) ./ sqrt(length(tmpcp(session.bins >=(0.05+addPulseTime(iPulse)) & session.bins <= (0.2+addPulseTime(iPulse)))));
%         
%         errorbar( ppk.m(iPulse), cpPulse.m(iPulse), ppk.se(iPulse), 'horizontal', 'o', 'Color', clr(iPulse, :), 'MarkerFaceColor', clr(iPulse, :), 'MarkerSize', 6);
%         errorbar( ppk.m(iPulse), cpPulse.m(iPulse), cpPulse.se(iPulse), 'vertical', 'o', 'Color', clr(iPulse, :), 'MarkerFaceColor', clr(iPulse, :), 'MarkerSize', 6);
%         
%     end % pulses
%     
%     legend('p1', '', 'p2', '', 'p3', '', 'p4', '', 'p5', '', 'p6', '', 'p7', '')
%     xlabel('Psychophysical Weight', 'FontSize', 12)
%     ylabel('Pulse CP', 'FontSize', 12)
%     xlim ([0.05 0.6]);
%     ylim ([0.48 0.57]);
%     axis('square')
%     supertitle('Direction model', 18)
%     
%     
%     % pta peaks vs. ppk
%     figure(2)
%     if cond==1
%         tmppta = max(flat.dir.combo.mPTA);
%         subplot(1,3,2); hold on
%         axis('square')
%         title('Flat', 'FontSize', 15)
%     end
%     if cond==2
%         tmppta = max(late.dir.combo.mPTA);
%         subplot(1,3,3); hold on
%         axis('square')
%         title('Late', 'FontSize', 15)
%     end
%     if cond==3
%         tmppta = max(early.dir.combo.mPTA);
%         subplot(1,3,1); hold on
%         axis('square')
%         title('Early', 'FontSize', 15)        
%     end
%     
%     for iPulse = 1:7
%         plot(ppk.m(iPulse), tmppta(iPulse), 'o', 'Color', clr(iPulse, :), 'MarkerFaceColor', clr(iPulse, :), 'MarkerSize', 6);
%     end    % pulses
%     
%     legend('p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7')
%     xlabel('Psychophysical Weight', 'FontSize', 12)
%     ylabel('PTA peak', 'FontSize', 12)
%     xlim([0.1 0.6]);
% end % three conditions
% %% choice model scatter plots
% 
% figure(3)
% figure(4)
% % pulse cp vs. ppk
% % get avg cp during a 150ms pulse
% for cond = 1:3
%     
%     figure(3)
%     if cond == 1
%         tmpcp  = flat.cho.cpm;        
%         ppk.m  = mean(flat.cho.combo.ppk);
%         ppk.se = std(flat.cho.combo.ppk) / sqrt(size(flat.cho.combo.ppk, 1));
%         
%         subplot(1,3,2); hold on
%         clr = winter(7);
%         title('Flat', 'FontSize', 15)
%     end
%     
%     if cond == 2
%         tmpcp = late.cho.cpm;
%         ppk.m = mean(late.cho.combo.ppk);
%         ppk.se = std(late.cho.combo.ppk) / sqrt(size(flat.cho.combo.ppk, 1));
%         
%         subplot(1,3,3); hold on
%         clr = autumn(7);
%         title('Late', 'FontSize', 15)
%     end
%     
%     if cond == 3
%         tmpcp = early.cho.cpm;
%         ppk.m = mean(early.cho.combo.ppk);
%         ppk.se = std(early.cho.combo.ppk) / sqrt(size(flat.cho.combo.ppk, 1));
%         
%         subplot(1,3,1); hold on
%         clr = hot(12);
%         title('Early', 'FontSize', 15)
%     end
%     
%     addPulseTime = [0 .150 .300 .450 .600 .750 .900];
%     for iPulse = 1:7
%         
%         cpPulse.m(iPulse) = mean(tmpcp(session.bins >=(0.05+addPulseTime(iPulse)) & session.bins <= (0.2+addPulseTime(iPulse))));
%         cpPulse.se(iPulse) = std(tmpcp(session.bins >=(0.05+addPulseTime(iPulse)) & session.bins <= (0.2+addPulseTime(iPulse)))) ./ sqrt(length(tmpcp(session.bins >=(0.05+addPulseTime(iPulse)) & session.bins <= (0.2+addPulseTime(iPulse)))));
%         
%         errorbar( ppk.m(iPulse), cpPulse.m(iPulse), ppk.se(iPulse), 'horizontal', 'o', 'Color', clr(iPulse, :), 'MarkerFaceColor', clr(iPulse, :), 'MarkerSize', 6);
%         errorbar( ppk.m(iPulse), cpPulse.m(iPulse), cpPulse.se(iPulse), 'vertical', 'o', 'Color', clr(iPulse, :), 'MarkerFaceColor', clr(iPulse, :), 'MarkerSize', 6);
%         
%     end % pulses
%     
%     legend('p1', '', 'p2', '', 'p3', '', 'p4', '', 'p5', '', 'p6', '', 'p7', '')
%     xlabel('Psychophysical Weight', 'FontSize', 12)
%     ylabel('Pulse CP', 'FontSize', 12)
%     xlim ([0.05 0.6]);
%     ylim ([0.53 0.595]);
%     axis('square')
%     supertitle('Choice model', 18)
%     
%     
%     % pta peaks vs. ppk
%     figure(4)
%     if cond==1
%         tmppta = max(flat.cho.combo.mPTA);
%         subplot(1,3,2); hold on
%         axis('square')
%         title('Flat', 'FontSize', 15)
%     end
%     if cond==2
%         tmppta = max(late.cho.combo.mPTA);
%         subplot(1,3,3); hold on
%         axis('square')
%         title('Late', 'FontSize', 15)
%     end
%     if cond==3
%         tmppta = max(early.cho.combo.mPTA);
%         subplot(1,3,1); hold on
%         axis('square')
%         title('Early', 'FontSize', 15)        
%     end
%     
%     for iPulse = 1:7
%         plot(ppk.m(iPulse), tmppta(iPulse), 'o', 'Color', clr(iPulse, :), 'MarkerFaceColor', clr(iPulse, :), 'MarkerSize', 6);
%     end    % pulses
%     
%     legend('p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7')
%     xlabel('Psychophysical Weight', 'FontSize', 12)
%     ylabel('PTA peak', 'FontSize', 12)
%     xlim([0.1 0.6]);
% end % three conditions
