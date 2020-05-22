% plot decoding by subject

leo.flat.dir  = decodingConditionAverage('direction', {'leo', 'flat'}, true);
leo.late.dir  = decodingConditionAverage('direction', {'leo', 'late'}, true);
leo.early.dir = decodingConditionAverage('direction', {'leo', 'early'}, true);

leo.flat.cho = decodingConditionAverage('choice', {'leo', 'flat'}, true);
leo.late.cho = decodingConditionAverage('choice', {'leo', 'late'}, true);
leo.early.cho = decodingConditionAverage('choice', {'leo', 'early'}, true);


nancy.flat.dir  = decodingConditionAverage('direction', {'nancy', 'flat'}, true);
nancy.late.dir  = decodingConditionAverage('direction', {'nancy', 'late'}, true);
nancy.early.dir = decodingConditionAverage('direction', {'nancy', 'early'}, true);

nancy.flat.cho = decodingConditionAverage('choice', {'nancy', 'flat'}, true);
nancy.late.cho = decodingConditionAverage('choice', {'nancy', 'late'}, true);
nancy.early.cho = decodingConditionAverage('choice', {'nancy', 'early'}, true);

%%
% cpt
%figure

figure
%subplot(234); hold on
subplot(231); hold on
leo.early.dir.cpm  = smooth( mean(leo.early.dir.cpt), 10);
leo.early.dir.cpse = smooth( (std(leo.early.dir.cpt) / sqrt(size(leo.early.cho.cpt, 1))), 9);
boundedline(session.bins, leo.early.dir.cpm, leo.early.dir.cpse, 'cmap', [.8 0 0], 'alpha');
%boundedline(session.bins_fromGo, dS.early.dir.cpm, dS.early.dir.cpse, 'cmap', [.8 0 0], 'alpha');

leo.early.cho.cpm  = smooth( mean(leo.early.cho.cpt), 10);
leo.early.cho.cpse = smooth( (std(leo.early.cho.cpt) / sqrt(size(leo.early.cho.cpt, 1))), 9);
boundedline(session.bins, leo.early.cho.cpm, leo.early.cho.cpse, 'cmap', [.3 0 0], 'alpha');
%boundedline(session.bins_fromGo, dS.early.cho.cpm, dS.early.cho.cpse, 'cmap', [.3 0 0], 'alpha');

axis square
ylabel('Choice probability')

%subplot(235); hold on
subplot(232); hold on
leo.flat.dir.cpm  = smooth( mean(leo.flat.dir.cpt), 10);
leo.flat.dir.cpse = smooth( (std(leo.flat.dir.cpt) / sqrt(size(leo.flat.cho.cpt, 1))), 9);
boundedline(session.bins, leo.flat.dir.cpm, leo.flat.dir.cpse, 'cmap', [0 0 .8], 'alpha');
%boundedline(session.bins_fromGo, dS.flat.dir.cpm, dS.flat.dir.cpse, 'cmap', [0 0 .8], 'alpha');

leo.flat.cho.cpm  = smooth( mean(leo.flat.cho.cpt), 10);
leo.flat.cho.cpse = smooth( (std(leo.flat.cho.cpt) / sqrt(size(leo.flat.cho.cpt, 1))), 9);
boundedline(session.bins, leo.flat.cho.cpm, leo.flat.cho.cpse, 'cmap', [0 0 .3], 'alpha');
%boundedline(session.bins_fromGo, dS.flat.cho.cpm, dS.flat.cho.cpse, 'cmap', [0 0 .3], 'alpha');

%xlim([-.5 1.5])
%ylim([0.46 .65])
axis square
xlabel('time (from GO)')

%subplot(236); hold on
subplot(233); hold on
leo.late.dir.cpm  = smooth( mean(leo.late.dir.cpt), 10);
leo.late.dir.cpse = smooth( (std(leo.late.dir.cpt) / sqrt(size(leo.late.cho.cpt, 1))), 9);
boundedline(session.bins, leo.late.dir.cpm, leo.late.dir.cpse, 'cmap', [.8 .8 0], 'alpha');
%boundedline(session.bins_fromGo, dS.late.dir.cpm, dS.late.dir.cpse, 'cmap', [.8 .8 0], 'alpha');

leo.late.cho.cpm  = smooth( mean(leo.late.cho.cpt), 10);
leo.late.cho.cpse = smooth( (std(leo.late.cho.cpt) / sqrt(size(leo.late.cho.cpt, 1))), 9);
boundedline(session.bins, leo.late.cho.cpm, leo.late.cho.cpse, 'cmap', [.4 .4 0], 'alpha');
%boundedline(session.bins_fromGo, dS.late.cho.cpm, dS.late.cho.cpse, 'cmap', [.4 .4 0], 'alpha');

%xlim([-.5 1.5])
%ylim([0.46 .65])
axis square


subplot(234)
h=histogram(0-leo.early.dir.timeToGo,10);
h.Normalization = 'probability';
axis square
xlim([-1 0]);
h.FaceColor= [.7 0 0];

subplot(2,3,5)
h=histogram(0-leo.flat.dir.timeToGo,10);
h.Normalization = 'probability';
axis square
xlim([-1 0]);
h.FaceColor= [0 0 .8];

subplot(2,3,6)
h=histogram(0-leo.late.dir.timeToGo,10);
h.Normalization = 'probability';
axis square
xlim([-1 0]);
h.FaceColor= [.7 .7 0];
%%

figure
%subplot(234); hold on
subplot(231); hold on
nancy.early.dir.cpm  = smooth( mean(nancy.early.dir.cpt), 10);
nancy.early.dir.cpse = smooth( (std(nancy.early.dir.cpt) / sqrt(size(nancy.early.cho.cpt, 1))), 9);
boundedline(session.bins, nancy.early.dir.cpm, nancy.early.dir.cpse, 'cmap', [.8 0 0], 'alpha');
%boundedline(session.bins_fromGo, dS.early.dir.cpm, dS.early.dir.cpse, 'cmap', [.8 0 0], 'alpha');

nancy.early.cho.cpm  = smooth( mean(nancy.early.cho.cpt), 10);
nancy.early.cho.cpse = smooth( (std(nancy.early.cho.cpt) / sqrt(size(nancy.early.cho.cpt, 1))), 9);
boundedline(session.bins, nancy.early.cho.cpm, nancy.early.cho.cpse, 'cmap', [.3 0 0], 'alpha');
%boundedline(session.bins_fromGo, dS.early.cho.cpm, dS.early.cho.cpse, 'cmap', [.3 0 0], 'alpha');

axis square
ylabel('Choice probability')

%subplot(235); hold on
subplot(232); hold on
nancy.flat.dir.cpm  = smooth( mean(nancy.flat.dir.cpt), 10);
nancy.flat.dir.cpse = smooth( (std(nancy.flat.dir.cpt) / sqrt(size(nancy.flat.cho.cpt, 1))), 9);
boundedline(session.bins, nancy.flat.dir.cpm, nancy.flat.dir.cpse, 'cmap', [0 0 .8], 'alpha');
%boundedline(session.bins_fromGo, dS.flat.dir.cpm, dS.flat.dir.cpse, 'cmap', [0 0 .8], 'alpha');

nancy.flat.cho.cpm  = smooth( mean(nancy.flat.cho.cpt), 10);
nancy.flat.cho.cpse = smooth( (std(nancy.flat.cho.cpt) / sqrt(size(nancy.flat.cho.cpt, 1))), 9);
boundedline(session.bins, nancy.flat.cho.cpm, nancy.flat.cho.cpse, 'cmap', [0 0 .3], 'alpha');
%boundedline(session.bins_fromGo, dS.flat.cho.cpm, dS.flat.cho.cpse, 'cmap', [0 0 .3], 'alpha');

%xlim([-.5 1.5])
%ylim([0.46 .65])
axis square
xlabel('time (from GO)')

%subplot(236); hold on
subplot(233); hold on
nancy.late.dir.cpm  = smooth( mean(nancy.late.dir.cpt), 10);
nancy.late.dir.cpse = smooth( (std(nancy.late.dir.cpt) / sqrt(size(nancy.late.cho.cpt, 1))), 9);
boundedline(session.bins, nancy.late.dir.cpm, nancy.late.dir.cpse, 'cmap', [.8 .8 0], 'alpha');
%boundedline(session.bins_fromGo, dS.late.dir.cpm, dS.late.dir.cpse, 'cmap', [.8 .8 0], 'alpha');

nancy.late.cho.cpm  = smooth( mean(nancy.late.cho.cpt), 10);
nancy.late.cho.cpse = smooth( (std(nancy.late.cho.cpt) / sqrt(size(nancy.late.cho.cpt, 1))), 9);
boundedline(session.bins, nancy.late.cho.cpm, nancy.late.cho.cpse, 'cmap', [.4 .4 0], 'alpha');
%boundedline(session.bins_fromGo, dS.late.cho.cpm, dS.late.cho.cpse, 'cmap', [.4 .4 0], 'alpha');

%xlim([-.5 1.5])
%ylim([0.46 .65])
axis square


subplot(234)
h=histogram(0-nancy.early.dir.timeToGo,10);
h.Normalization = 'probability';
axis square
xlim([-1 0]);
h.FaceColor= [.7 0 0];
ylabel('proportion trials')

subplot(2,3,5)
h=histogram(0-nancy.flat.dir.timeToGo,10);
h.Normalization = 'probability';
axis square
xlim([-1 0]);
h.FaceColor= [0 0 .8];
xlabel('stimulus off time (from GO)')

subplot(2,3,6)
h=histogram(0-nancy.late.dir.timeToGo,10);
h.Normalization = 'probability';
axis square
xlim([-1 0]);
h.FaceColor= [.7 .7 0];