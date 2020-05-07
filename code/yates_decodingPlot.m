addpath analysis_scripts\decoding\
figDir='figures/decodingMTredux2019';
dataPath = getpref('mtlipglm', 'dataPath');

decoder = 'logistic-L2-glmnet';
% decoder = 'linear-lda';
% decoder = 'logistic-Elastic-Net-glmnet';

directory = 'data/decoding_redux';

fid = fopen(fullfile(figDir, sprintf('stats2_%s.txt', decoder)), 'w');

%% Count decoding

[experiments,D] = getInclusionCriterion(5, fid);

%% load up count decoding
date = ['_20180608'];
% date = '';
clear Sblind Nblind Sopt nOpt Squad Nquad
for kEx = 1:numel(D)
    
    fname = fullfile(directory, ...
        sprintf('Count_%s_%s%s.mat', D(kEx).exname, decoder, date));
    
    tmp = load(fname);
    
    Sopt(kEx) = tmp.Sopt;
    Nopt(kEx) = tmp.Nopt;
    Sblind(kEx) = tmp.Sblind;
    Nblind(kEx) = tmp.Nblind;
    if isfield(tmp, 'Squad')
        Squad(kEx) = tmp.Squad;
        Nquad(kEx) = tmp.Nquad;
    end
    
end



%% print out statistics 
z_units = [];
for kEx = 1:numel(D)

ix = Sopt(kEx).trialIdx;
x = D(kEx).stim.dir(ix);     % direction
y = D(kEx).stim.choice(ix);  % choice
c = x==y;                    % correct
n = numel(y);                % number of trials

% bootstrap
B = 1e3; % number
bi = randi(n, [n B]); % indices
pcb = mean(c(bi));    % boot percent correct
mu = mean(pcb);
sd = std(pcb);

z_units = [z_units (Sopt(kEx).pcIndividualCovariates - mu) / sd]; %#ok<AGROW>
z_acc(kEx) = (Sopt(kEx).decodingAccuracy - mu) / sd;
end


fprintf(fid, 'Compare absolute performance of single units to monkey at p<.05 (two-sided)\n');
fprintf(fid, '%d/%d units better than monkey. %d/%d units worse.\n', sum(z_units > 1.96), numel(z_units), sum(z_units < -1.96), numel(z_units));
fprintf(fid, '%d/%d sessions better than monkey. %d/%d sessions worse.\n', sum(z_acc > 1.96), numel(z_acc), sum(z_acc < -1.96), numel(z_acc));

%% example session nmf pmf
% figure(1); clf
% kEx = strcmp({D.exname}, 'n20150324a');
% plot(Nopt(kEx).pmf.x, Nopt(kEx).pmf.y, 'ok', 'MarkerFaceColor', 'k'); hold on
% h = plot(Nopt(kEx).pmf.xx, Nopt(kEx).pmf.yy, '-k');
% 
% nid = max(Nopt(kEx).nmf.dprime)==Nopt(kEx).nmf.dprime;
% clr = repmat(.5, 1, 3);
% plot(Nopt(kEx).nmf.x, Nopt(kEx).nmf.y(:,nid), 'o', 'Color', clr, 'MarkerFaceColor', clr); hold on
% h(2) = plot(Nopt(kEx).nmf.xx(:,nid), Nopt(kEx).nmf.yy(:,nid), '-', 'Color', clr);
% 
% plot(Nopt(kEx).popnmf.x, Nopt(kEx).popnmf.y, 'or', 'MarkerFaceColor', 'r'); hold on
% h(3) = plot(Nopt(kEx).popnmf.xx, Nopt(kEx).popnmf.yy, '-r');
% set(gca, 'YTick', .5:.1:1)
% pdsa.offsetAxes(gca)
% 
% pdsa.fixfigure(gcf, 8, [2 2])
% 
% legend(h, {'Monkey', 'Best Neuron', 'Population'}, 'Location', 'SouthEast', 'Box', 'off', 'FontSize', 8, 'FontName', 'Arial')
% xlabel('Motion Strength')
% ylabel('Proportion Correct')
% saveas(gcf, fullfile(figDir, 'fig04_panelA.pdf'))
N = Nopt;
exampleSessions={'p20140304', 'p20140305', 'n20150316c', 'n20150324a'};
nExamples = numel(exampleSessions);
sx = ceil(sqrt(nExamples));
sy = round(sqrt(nExamples));

fig3 = figure(128); clf
clr = 'r';
for k = 1:nExamples
    
    subplot(sy,sx,k)
    
    kEx = find(strcmp(exampleSessions{k}, experiments));
    D(kEx).plotPMF(N(kEx).popnmf, 'r')
    D(kEx).plotPMF(N(kEx).pmf, 'k')
    [~,id]=min(N(kEx).nmf.alpha);
    plot(N(kEx).nmf.x, N(kEx).nmf.y(:,id), '.', 'Color', .5*[1 1 1])
    plot(N(kEx).nmf.xx(:,id), N(kEx).nmf.yy(:,id), 'Color', .5*[1 1 1])
    ylim([.4 1])
    xlabel('Motion Strength')
    ylabel('Proportion Correct')
    xlim([0 1.65])
    
    if kEx==1
        set(gca, 'Xtick', 0:.5:1.5)
        ylim([.5 1])
        text(.6,.65, 'Monkey')
        text(.6,.6, 'Population')
        text(.6,.55, 'Best Neuron')
        plot(.5+[0 .05], .65*[1 1], 'k')
        plot(.5+[0 .05], .6*[1 1], 'Color', clr)
        plot(.5+[0 .05], .55*[1 1], 'Color', .5*[1 1 1])
    end
    
    title('')
end

pdsa.fixfigure(fig3, 8, [5.5 5.5])
saveas(fig3, fullfile(figDir, 'fig03_examples-full-ver01.pdf'))
%% plot count decoding

% --- NmfPmfSummary contains the analyses and figures required to produce
% figure 4 (except for panel A, which is produced above)
fig4 = figure_NmfPmfSummary(D, Nopt, Sopt,1, fid);

pdsa.fixfigure(fig4, 8, [7.5 3], 'FontName', 'Arial')
saveas(fig4, fullfile(figDir, sprintf('fig04_nmf-count-ver01-%s.pdf', decoder)))

fig4 = figure_NmfPmfSummary(D, Nopt, Sopt,0, fid);

pdsa.fixfigure(fig4, 8, [7.5 2.5], 'FontName', 'Arial')
saveas(fig4, fullfile(figDir, sprintf('fig04_nmf-count-ver02-%s.pdf', decoder)))

fig4 = figure_NmfPmfSummary(D, Nopt, Sopt,2, fid);

pdsa.fixfigure(fig4, 8, [5.5 5.5], 'FontName', 'Arial')
saveas(fig4, fullfile(figDir, sprintf('fig04_nmf-count-ver03-%s.pdf', decoder)))

% --- CompareDecoders contains the analyses required to produce figure 3
fig5 = figure_CompareDecoders(D, Sopt, Sblind,fid);

pdsa.fixfigure(fig5(1), 8, [5.5 5.5], 'FontName', 'Arial')
saveas(fig5(1), fullfile(figDir, sprintf('fig03_count-correlations-ver03-%s.pdf', decoder)))

% regenerate entire figure just to get the quadratic decoder plot for panel D
fig6 = figure_CompareDecoders(D, Sopt, Squad,fid);

pdsa.fixfigure(fig6, 8, [5.5 5.5], 'FontName', 'Arial')
saveas(fig6, fullfile(figDir, sprintf('fig03_panelD_opt_quad_compare-%s.pdf', decoder)))

%%

fprintf(fid, '------------------------------------------------------------\n');
fprintf(fid, '------------------------------------------------------------\n');
fprintf(fid, '------------------------------------------------------------\n\n\n');
fprintf(fid, 'Instantaneous Decoder Analyeses\n\n\n\n');
%% Instantaneous decoder
zscoreRegressors = true;

clear D
for kEx = 1:numel(experiments)
    exname = experiments{kEx};
    
    stim = getStim(exname);
    neurons = getNeurons(exname, dataPath);
    neurons(~[neurons(:).isMT]) = []; % only MT
    
    D(kEx) = decode_session(stim, neurons, ...
        'winPre', -.3, ...
        'winPost', .4); %#ok<*SAGROW>
end

%% plot psth
ctr = 0;

n = min(arrayfun(@(x) numel(x.spikeBinCenters), D));
nC = 5;

mu = nan(n, nC, sum(arrayfun(@(x) x.nFitNeurons, D)));

for kEx = 1:numel(D)
    [~, ind] = sort(D(kEx).stim.coh);
    for kNeuron = 1:D(kEx).nFitNeurons
        spks = D(kEx).spikes(ind,:,kNeuron);
        
        if D(kEx).dPrimeNeuron(kNeuron)<0
            spks = flipud(spks);
        end
        
        bc = ones(5,1)/5;
        smspks = filter(bc, 1, spks')';
        
        nT = size(spks, 1);
        
        binedges = round(linspace(1,nT, nC+1));
        
        
        for i = 1:nC
            mu(:,i,kNeuron+ctr) = mean(smspks(binedges(i):binedges(i+1),1:n));
        end
        
    end
    ctr = ctr + D(kEx).nFitNeurons;
end
figure(1); clf
t = D(1).spikeBinCenters(1:n);
cmap = pdsa.cbrewer('jake', 'rdbu', nC);
set(gcf, 'DefaultAxesColorOrder', cmap)
plot(t, mean(mu,3)/mean(diff(t)))
xlim([-.2 1.5])
xlabel('Time from motion onset (s)')
ylabel('Firing rate')
text(1, 35, sprintf('n=%d', size(mu,3)))
pdsa.fixfigure(gcf, 8, [2.5 2.5], 'FontName', 'Arial')
saveas(gcf, fullfile(figDir, 'fig01_poppsth.pdf'))

%%



%% load up decoding results
decoder = 'logistic-L2-glmnet';

clear Scount Sinst Scum dcComp dcFull
Scum  = [];
Sinst = [];
dcComp = [];
dcFull = [];
nBins = 173;

for kEx = 1:numel(D)
    
    fname = fullfile(directory, ...
        sprintf('Inst_%s_%s.mat', D(kEx).exname, decoder));
    
    tmp = load(fname);
    
    dcComp = [dcComp; tmp.dcComp(1:nBins)];
    Scum = [Scum; tmp.Scum(1:nBins)];
%     if isfield(tmp.Sinst, 'yhatBoot')
%         for i = 1:numel(tmp.Sinst)
%             tmp.Sinst(i).yhat = tmp.Sinst(i).yhatBoot;
%         end
%     end
    if isfield(tmp.Sinst, 'yhatBoot')
    tmp.Sinst = rmfield(tmp.Sinst, 'yhatBoot');
    tmp.Sinst = rmfield(tmp.Sinst, 'ycorr');
    end
    Sinst = [Sinst; tmp.Sinst(1:nBins)];
    dcFull = [dcFull; tmp.dcFull];
end

%% plotting

pcFixed = cell2mat(arrayfun(@(x) x.accuracyFull(1:nBins), dcFull, 'UniformOutput', false))';
pcInst  = cell2mat(arrayfun(@(x) x.accuracyInst2(1:nBins), dcFull, 'UniformOutput', false))';


fig = figure(122); clf

% --------------------------------------------------------------------------
% --- Instantaneous decoder vs. Fixed over time
ax = subplot('Position', [0.100 0.16 0.2134 0.76]);

xat = D(1).spikeBinCenters(1:nBins);
cmap = lines;

m = mean(pcInst,2);
s = std(pcInst, [], 2)/sqrt(numel(D));
clr = cmap(3,:);
errorbarFill(xat, m, s, 'k', 'FaceColor', clr, 'FaceAlpha', .5, 'EdgeColor', 'none'); hold on
plot(xat, m, 'Color', clr); hold on
% plot(xat, m+s, '--', 'Color', clr)
% plot(xat, m-s, '--', 'Color', clr)
xlabel('Time')
ylabel('Accuracy')

plot(.4, .74, '.', 'Color', clr);
text(.5, .74, 'Instantaneous')

m = mean(pcFixed,2);
s = std(pcFixed, [], 2)/sqrt(numel(dcFull));
clr = cmap(4,:);
errorbarFill(xat, m, s, 'k', 'FaceColor', clr, 'FaceAlpha', .5, 'EdgeColor', 'none'); hold on
plot(xat, m, 'Color', clr); hold on
% plot(xat, m+s, '--', 'Color', clr)
% plot(xat, m-s, '--', 'Color', clr)

plot(.4, .72, '.', 'Color', clr);
text(.5, .72, 'Fixed')

xlim([-.2 1.3])

set(ax, 'Box', 'off')


binidx = xat>0 & xat < 1.2;
rat = pcFixed(binidx,:) ./ pcInst(binidx,:);

gm = geomean(rat(randi(numel(rat), 2e3)));

fprintf(fid, 'Geometric Mean: %02.4f +- %02.4f (bootstrapped CI) \n', mean(gm), 2*std(gm));

mu0 = 1;
[h,pval,~,tstats] = ttest(rat(:), mu0);
if h
    fprintf(fid, 'Reject the null hypothesis that Instantaneous and Fixed decoders are the same\n');
else
    fprintf(fid, 'Cannot reject the null hypothesis that Instantaneous and Fixed are the same\n');
end
fprintf(fid, 'two-sided ttest. t(%d) = %02.4f, p = %d\n', tstats.df, tstats.tstat, pval);


[pval, hyp, stats] = signrank(rat(:), mu0);
if hyp
    fprintf(fid, 'Reject the Null\n');
else
    fprintf(fid, 'Cannot reject the null\n');
end
    
fprintf(fid, 'two-sided Wilcoxin signed rank test. signrank(%d) = %02.2f, p = %d\n', numel(rat)-1, stats.signedrank, pval);
% --------------------------------------------------------------------------
% --- stabilitiy of decoding
ax = subplot('Position', [0.4 0.1600 0.3 0.76]);

% The stability of category
pcl = [];
for t = 1:nBins
    pcLagged = cell2mat(arrayfun(@(x) x.accuracyLagged(1:nBins), dcComp(:,t), 'UniformOutput', false));
    
    pcl = [pcl; mean(pcLagged)];
end
    

imagesc(xat, xat, pcl); colormap jet
axis xy
h = colorbar;
h.Label.String = '% correct';

xlabel('Time Tested')
ylabel('Time Trained')
set(ax, 'Box', 'off', 'XTick', 0:.5:1, 'YTick', 0:.5:1)

% --------------------------------------------------------------------------
% --- weights over time
ax = subplot('Position', [0.79 0.16 0.1819 0.76]);

S = Sinst(:,1:nBins);

w  = [];
pd = [];
for kEx = 1:size(S,1)
    wtmp = cell2mat(arrayfun(@(x) x.wts', S(kEx,:)', 'UniformOutput', false));
    w  = [w wtmp];
    pd = [pd D(kEx).prefDir(:)'];
end

w(:,isnan(pd)) = [];
pd(isnan(pd)) = [];

bs = 30;
binEdge = -180:bs:180;
[~, id] = histc(pd, binEdge);

n = numel(unique(id));
I = zeros(n, nBins);
for i = 1:n
   iix = id==i;
   I(i,:) = mean(w(:,iix),2)';
end


binCenter = binEdge(1:end-1)+bs/2;
imagesc(xat(1:nBins), binCenter, I, 2*[-1 1])
axis xy
set(ax, 'YTick', -180:90:180)
colormap jet
xlabel('Time')
ylabel('Preferred Direction')


pdsa.fixfigure(fig, 8, [7.5 2])
saveas(fig, fullfile(figDir, 'fig06_inst-decoder-ver01.pdf'))

% colormap(gray.^1.5)
pdsa.fixfigure(fig, 8, [7.5 2])
saveas(fig, fullfile(figDir, 'fig06_inst-decoder-ver02.pdf'))

%% compare first half to second half decoding

%% load up count decoding
clear transient sustained
for kEx = 1:numel(D)
    
    fname = fullfile(directory, ...
        sprintf('EarlyLate2_%s_%s.mat', D(kEx).exname, decoder));
    
    tmp = load(fname);
    
    transient(kEx) = tmp.early;
    sustained(kEx) = tmp.late;
end


%% Peak sensitivity of temporal decoder


[~, peakid]=max(pcInst);

fig7 = figure(123); clf

% --------------------------------------------------------------------------
% Time until peak
ax = subplot('Position', [0.1300 0.15 0.2134 0.76]);
histogram(xat(peakid), 'binEdges', 0:.05:1.2, 'Normalization', 'probability', 'FaceColor', .5*[1 1 1], 'EdgeColor', 'k')
xlabel('Time until peak accuracy')
ylabel('Fraction of Sessions')

xlim([0 1])
set(ax, 'Box', 'off', 'XTick', 0:.25:1, 'YTick', 0:.1:.5)


% --------------------------------------------------------------------------
% --- compare first half to second half decoding
ax = subplot('Position', [0.4108 0.15 0.2134 0.76]);
x   = [sustained(:).decodingAccuracy];
xl  = arrayfun(@(x) x.decodingAccuracyCi(1), sustained(:))';
xu  = arrayfun(@(x) x.decodingAccuracyCi(2), sustained(:))';
y   = [transient(:).decodingAccuracy];
yl  = arrayfun(@(x) x.decodingAccuracyCi(1), transient(:))';
yu  = arrayfun(@(x) x.decodingAccuracyCi(2), transient(:))';

plot([x; x], [y + yl; y + yu], 'Color', .5*[1 1 1]); hold on
plot([x+xl; x+xu], [y; y], 'Color', .5*[1 1 1]);

cmap=hot(100);
cmap=cmap(1:65,:);
colormap(cmap)

n = arrayfun(@(x) sum(x.neuronIx), D);
plotWeights([x(:),y(:)],n(:), 10);
axis on
axis normal
% plot(x,y, 'k.')
hold on
xd = [.6 1];
plot(xd, xd, 'k')
xlim(xd); ylim(xd)
ylabel('Accuracy (First Half)')
xlabel('Accuracy (Second Half)')

set(ax, 'Box', 'off', 'XTick', .6:.1:1, 'YTick', .6:.1:1)

% --------------------------------------------------------------------------
% --- Psychophysical kernel
ax = subplot('Position', [0.6916 0.15 0.2134 0.76]);

X = cell2mat(arrayfun(@(x) x.stim.pulses(x.stim.frozenTrialIx==0,:), D(:), 'uniformOutput', false));
y = cell2mat(arrayfun(@(x) x.stim.choice(x.stim.frozenTrialIx==0)==1, D(:), 'uniformOutput', false));


[b, ~, stats] = glmfit(X,y, 'binomial');
b = b(2:end);
ptime = 1:7;
errorbarFill(ptime, b, stats.se(2:end), 'k', 'FaceColor', 'k', 'FaceAlpha', .5, 'EdgeColor', 'none'); hold on
plot(ptime, b, 'k-')
hold on
plot(ptime, b, 'k.', 'MarkerSize', 10);
% plot(ptime, b+stats.se(2:end), 'k--')
% plot(ptime, b-stats.se(2:end), 'k--')


X = cell2mat(arrayfun(@(x,y) x.stim.pulses(y.trialIdx,:), D(:), Sopt(:), 'UniformOutput', false));
y = cell2mat(arrayfun(@(x) x.yhat==1, Sopt(:), 'UniformOutput', false));

[b, ~, stats] = glmfit(X,y, 'binomial');
b = b(2:end);
ptime = 1:7;

errorbarFill(ptime, b, stats.se(2:end), 'b', 'FaceColor', 'b', 'FaceAlpha', .5, 'EdgeColor', 'none'); hold on

plot(ptime, b, 'b-')
hold on

plot(ptime, b, 'b.', 'MarkerSize', 10);
% plot(ptime, b+stats.se(2:end), 'b--')
% plot(ptime, b-stats.se(2:end), 'b--')

xlabel('Pulse #')
ylabel('Neurometric Weight')
xlim([.5 7.5])
ylim([0 .8])
set(ax, 'Box', 'off', 'XTick', ptime, 'YTick', 0:.2:.8)
text(5, .9*max(ylim), sprintf('n = %d', numel(y)))

pdsa.fixfigure(fig7, 8, [7.5 2])
saveas(fig7, fullfile(figDir, 'fig07_early-sensitivity-ver01.pdf'))

%% --------------------------------------------------------------------------
% --- Psychophysical kernel separated by monkey
figure(111); clf

monkeyIx = arrayfun(@(x) strcmp(x.exname(1), 'n'), D);

% sessions that the neurometric performance was comparable to monkey
gIx = arrayfun(@(x) x.popnmf.pc, Nopt)./arrayfun(@(x) x.pmf.pc, Nopt) > .95;
mIx = {monkeyIx==0, monkeyIx==1, monkeyIx==1 & gIx};
labels = {'P', 'N','N(good)'};

n = 2;
for i = 1:n
    subplot(1,n,i)
    iix = mIx{i};
    lbl = labels{i};
    X = cell2mat(arrayfun(@(x) x.stim.pulses(x.stim.frozenTrialIx==0,:), D(iix)', 'uniformOutput', false));
    y = cell2mat(arrayfun(@(x) x.stim.choice(x.stim.frozenTrialIx==0)==1, D(iix)', 'uniformOutput', false));
    
    
    [b, ~, stats] = glmfit(X,y, 'binomial');
    b = b(2:end);
    ptime = 1:7;
    plot(ptime, b, 'k-')
    hold on
    plot(ptime, b, 'k.', 'MarkerSize', 10);
    plot(ptime, b+stats.se(2:end), 'k--')
    plot(ptime, b-stats.se(2:end), 'k--')
    
    
    X = cell2mat(arrayfun(@(x,y) x.stim.pulses(y.trialIdx,:), D(iix)', Sopt(iix)', 'UniformOutput', false));
    y = cell2mat(arrayfun(@(x) x.yhat==1, Sopt(iix)', 'UniformOutput', false));
    
    [b, ~, stats] = glmfit(X,y, 'binomial');
    b = b(2:end);
    ptime = 1:7;
    plot(ptime, b, 'b-')
    hold on
    plot(ptime, b, 'b.', 'MarkerSize', 10);
    plot(ptime, b+stats.se(2:end), 'b--')
    plot(ptime, b-stats.se(2:end), 'b--')
    
    xlabel('Pulse #')
    ylabel('Weight')
    title(sprintf('Monkey %s', lbl))
    xlim([.5 7.5])
    % ylim([0 .8])
    set(gca, 'Box', 'off', 'XTick', ptime, 'YTick', 0:.2:max(ylim))
    
    text(6, .9*max(ylim), sprintf('n = %d', numel(y)))
    
    
end


pdsa.fixfigure(fig7, 8, [7.5 2])
saveas(fig7, fullfile(figDir, 'fig07b_PKmonkey.pdf'))




%% Load full temporal model
decoder = 'logistic-L2-glmnet';
% decoder = 'logistic-Elastic-Net-glmnet';

clear O Stmp Sfull Sb N Schoice
for kEx = 1:numel(experiments)
    
    exname = experiments{kEx};
    
    stim = getStim(exname);
    neurons = getNeurons(exname,dataPath);
    neurons(~[neurons(:).isMT]) = []; % only MT
    
    O(kEx) = decode_session(stim, neurons, 'binSize', .05); %#ok<*SAGROW>
    
    
    
    fname = fullfile(directory, ...
        sprintf('Full_%s_%s.mat', D(kEx).exname, decoder));
    
    tmp = load(fname);
    
    Stmp(kEx) = tmp.Stmp;
    Sfull(kEx) = tmp.Sfull;
    Sb(kEx) = tmp.Sb;
    
    fprintf(fid, '%s: %02.3f, %02.3f, %02.3f\n', O(kEx).exname, Stmp(kEx).decodingAccuracy, Sfull(kEx).decodingAccuracy, O(kEx).monkeyPC);
    
    fname = fullfile(directory, ...
        sprintf('FullChoice_%s_%s.mat', D(kEx).exname, decoder));
    
    tmp = load(fname);
    Schoice(kEx) = tmp.Sfull;
end

%% Plot cumulative performance
fig8 = figure(124); clf

nBins = 173;
pccum = cell2mat(arrayfun(@(x) x.decodingAccuracy, Scum(:,1:nBins), 'UniformOutput', false));
pccumCI = cell2mat(arrayfun(@(x) x.decodingAccuracyCI, Scum(:,1:nBins), 'UniformOutput', false));
xat = D(1).spikeBinCenters(1:nBins);

% -------------------------------------------------------------------------
% --- Cumulative performance

cmap=hot(100);
cmap=cmap(1:65,:);
colormap(cmap)

clrs = getColors(arrayfun(@(x) sum(x.neuronIx), D));
        
ax = subplot(1,2,1);
m = mean(pccum(2:end,:));
s = std(pccum);
% cmap = 1-lines*.88;
for i = 1:numel(D)
plot(xat, pccum(i,:), 'Color', clrs(i,:)); hold on
end
plot(xat(1:nBins), m, 'k', 'Linewidth', 2);
xlim([-.2 1.3])
xlabel('Time')
ylabel('Accuracy')

% calculate saturation time for each session
finalPerformance = pccum(:,end)-pccumCI(:,end)/2; % within SE of performance
nmtime = bsxfun(@rdivide, pccum', finalPerformance');

saturationTime = [];
for i = 1:numel(D)
    tcross = find(nmtime(:,i)>1, 1, 'first');
    if isempty(tcross)
        continue
    end
    saturationTime = [saturationTime; xat(tcross)];
    plot(xat(tcross), pccum(i,tcross), 'k.'); %, 'Color', clrs(i,:))
end

% plot histogram of cross times
ax2 = axes('Position', ax.Position);
bins = 0:.1:1.1;
cnt = histc(saturationTime,bins);
bar(bins, cnt, 'FaceColor', .5*ones(1,3))
ylim([0 10])
xlim([-.2 1.3])
set(ax2, 'box', 'off', 'YAxisLocation', 'right', 'Color', 'none')


ax = subplot(1,2,2);
% finalPerformance = pccum(:,end)-pccumCI(:,end)/2;
% 
% nmtime = bsxfun(@rdivide, pccum', finalPerformance');
nmtime = bsxfun(@rdivide, pccum', [D.monkeyPC]);
for i = 1:numel(D)
plot(xat, nmtime(:, i), 'Color', clrs(i,:)); hold on
end
plot(xat, mean(nmtime, 2), 'k', 'Linewidth', 2);
saturationTimeNP = [];
for i = 1:numel(D)
    tcross = find(nmtime(:,i)>1, 1, 'first');
    if isempty(tcross)
        continue
    end
    plot(xat(tcross), 1, 'ok'); hold on
    saturationTimeNP = [saturationTimeNP; xat(tcross)];
end
xlim([-.1 1.1])
plot(xlim, [1 1], 'k--')
ax2 = axes('Position', ax.Position);
bins = 0:.1:1.1;
cnt = histc(saturationTimeNP,bins);
xlim([-.1 1.1])

bar(ax2, bins, cnt, 'FaceColor', .5*[1 1 1])
set(ax2, 'box', 'off', 'YAxisLocation', 'right', 'Color', 'none')
% axis(ax2, 'off')
ylim(ax2, [0 20])
xlim(ax2, ax.XLim)
ylabel(ax2, '# sessions')

% ax = subplot(1,3,3);
% nmtime = bsxfun(@rdivide, pccum', [D.monkeyPC]);
% plot(xat, nmtime); hold on
% saturationTimeNP = [];
% for i = 1:numel(D)
%     tcross = find(nmtime(:,i)>1, 1, 'first');
%     if isempty(tcross)
%         continue
%     end
%     plot(xat(tcross), 1, 'ok'); hold on
%     saturationTimeNP = [saturationTimeNP; xat(tcross)];
% end
% 
% plot(xlim, [1 1], 'k--')
% ax2 = axes('Position', ax.Position);
% bins =  0:.1:1.2;
% cnt = histc(saturationTimeNP,bins);
% 
% bar(ax2, bins, cnt, 'FaceColor', .5*[1 1 1])
% % axis(ax2, 'off')
% ax2.Visible = 'off';
% ylim(ax2, [0 20])
% xlim(ax2, ax.XLim)

ylabel(ax, 'Neurometric/Psychometric')

pdsa.fixfigure(fig8, 8, [6 2])
saveas(fig8, fullfile(figDir, 'fig08a_cumulative.pdf'))
%% Temporal weights
% -------------------------------------------------------------------------
% --- Temporal weights strength
figure(2); clf
n = 20;
wT = [];
wC = [];

for kEx = 1:numel(D)
% nNeurons = sum(Stmp(kEx).wts~=0);
nTimeBins = numel(Sfull(kEx).binIdx);
W = reshape(Sfull(kEx).wts, nTimeBins, []);
Wcho = reshape(Schoice(kEx).wts, nTimeBins, []);



wTime = W*(mean(W))';
wChoTime = Wcho*(mean(Wcho))';

wT = [wT wTime(1:n)];
wC = [wC wChoTime(1:n)];

plot(O(kEx).spikeBinCenters(Sfull(kEx).binIdx), wTime, 'b');
hold on
plot(O(kEx).spikeBinCenters(Sfull(kEx).binIdx), wChoTime, 'r')

end

figure(3); clf

t = O(kEx).spikeBinCenters(1:n);
m = mean(wT,2);
s = std(wT, [],2)/sqrt(numel(D));
plot(t, m, 'b'); hold on
plot(t, m+s, 'b--')
plot(t, m-s, 'b--')

t = O(kEx).spikeBinCenters(1:n);
m = mean(wC,2);
s = std(wC, [],2)/sqrt(numel(D));
plot(t, m, 'r'); hold on
plot(t, m+s, 'r--')
plot(t, m-s, 'r--')

xlabel('Time (s)')
ylabel('Weight')

pdsa.fixfigure(gcf, 8, [3 2])
saveas(gcf, fullfile(figDir, 'fig9_temporal_weights.pdf'))

%%

return
%% Do we dare explore choice

figure(1); clf
for kEx = 1:numel(D)
    pcdir = mean(Sopt(kEx).yhat == D(kEx).stim.dir(Sopt(kEx).trialIdx));
    pccho = mean(Sopt(kEx).yhat == D(kEx).stim.choice(Sopt(kEx).trialIdx));
    
    plot(pcdir, pccho, '.'); hold on
end
plot(xd, xd, 'k')
xlabel('Accuracy (Direction Decoding)')
ylabel('Accuracy (Choice Decoding)')



1./(1+ -D(kEx).spikeCount(:,D(kEx).neuronIx)*Sopt(kEx).wts + Sopt(kEx).w0)

%% Continuoation of temporal weights

figure(3); clf
subplot(1,3,1)
plot(W)
subplot(1,3,2)
plot(Wcho)

subplot(1,3,3)
plot(O(kEx).spikeBinCenters(Sfull(kEx).binIdx), wTime, 'b');
hold on
plot(O(kEx).spikeBinCenters(Sfull(kEx).binIdx), wChoTime, 'r')

[u,s,v] = svd(W);

% plot(u(:,1)*sign(sum(u(:,1))))

figure(1); clf
subplot(1,2,1)
plot(W)
subplot(1,2,2)

sd = 3;
plot(u(:,1:sd)*s(1:sd,1:sd)*v(:,1:sd)')

plot(u(:,1:sd)*s(1:sd,1:sd)*sign(sum(u(:,1:sd)))')
hold on
plot(mean(bsxfun(@times, W, sign(sum(W))),2), 'k')


figure(2); clf
subplot(1,2,1)
imagesc(W)
subplot(1,2,2)
imagesc((bsxfun(@times, W, sign(sum(W)))))



%%
clf
w  = [];
pd = [];
nBins = min(arrayfun(@(x) numel(x.binIdx), Sfull));
for kEx = 1:numel(Stmp)
    neuronix = Stmp(kEx).wts~=0;
    pd_ = D(kEx).prefDir(neuronix);
    pd = [pd pd_(:)'];
    w_ = reshape(Sfull(kEx).wts, [], sum(neuronix));
    w  = [w w_(1:nBins,:)];
end

w(:,isnan(pd)) = [];
pd(isnan(pd)) = [];
[~, id] = sort(pd);


% wtmp = bsxfun(@rdivide, w, sum(w));
wtmp = bsxfun(@times, w, sign(sum(w)));
xat = D(1).spikeBinCenters(Sfull(1).binIdx);

m = mean(wtmp,2)';
s = std(wtmp, [],2)'/sqrt(size(wtmp,2));
xx = xat(1:nBins)-mean(diff(xat))/2;

plot([xx; xx], [m-s; m+s], '-', 'Color', .5*[1 1 1]); hold on
plot(xx, m, 'k.');
plot(xx, m, 'k-');

axis tight
plot(xlim, median(m)*[1 1 ], 'k--')
xlabel('Time (s)')
ylabel('|w|')
% stairs(xx, m+s, '-');
%%
% -------------------------------------------------------------------------
% --- Neurometric weights
figure(1); clf
patIx = cellfun(@(x) strcmp(x(1), 'p'), {O.exname});
patIx = true(size(patIx));
% patIx = 11;
X = cell2mat(arrayfun(@(x,y) x.stim.pulses(y.trialIdx,:), O(patIx)', Sfull(patIx)', 'UniformOutput', false));
y = cell2mat(arrayfun(@(x) x.yhat, Stmp(patIx)', 'UniformOutput', false));

qf = qfsmooth1D(size(X,2));
bb = decodeLogisticQuad(X,y, qf);

b = bb.wts;
se = bb.se;

ptime = 1:7;
ax = subplot(1,3,3);
plot(ptime, b, 'k-')
hold on
plot(ptime, b, 'k.', 'MarkerSize', 10)
plot(ptime, b+se, 'k--')
plot(ptime, b-se, 'k--')
xlabel('Pulse #')
ylabel('Neurometric Weight')
xlim([.5 7.5])
set(ax, 'Box', 'off', 'XTick', ptime, 'YTick', 0:.2:.8)

% psychophysical weights
y = cell2mat(arrayfun(@(x,y) x.stim.choice(y.trialIdx,:)==1, O(patIx)', Stmp(patIx)', 'UniformOutput', false));

bb = decodeLogisticQuad(X,y, 'qf', qf);
b = bb.wts;
se = bb.se;

plot(ptime, b, 'r-')
hold on
plot(ptime, b, 'r.', 'MarkerSize', 10)
plot(ptime, b+se, 'r--')
plot(ptime, b-se, 'r--')
xlabel('Pulse #')
ylabel('Neurometric Weight')
xlim([.5 7.5])
set(ax, 'Box', 'off', 'XTick', ptime, 'YTick', 0:.2:.8)

pdsa.fixfigure(fig8, 8, [7.5 2])
saveas(fig8, fullfile(figDir, 'fig08_optimal-temporalwts-ver01.pdf'))

%%
fclose(fid);

%%

withinError = pccum;
id = nan(size(withinError,1),1);
for i = 1:size(withinError,1)
    id(i) = find(withinError(i,:) > (pccum(i,end)-pccumCI(i,end)), 1);
    
%     figure(10); clf
%     plot(xat, withinError(i,:)); hold on
%     plot(xlim, pccum(i,end)*[1 1], 'r--')
%     plot(xat(id(i))*[1 1], ylim, 'k--')
%     pause
end

% th = .95;
% fracPerf = bsxfun(@rdivide, pccum', pccum(:,end)');
% id = nan(size(fracPerf,2),1);
% for i = 1:size(fracPerf,2)
%     id(i) = find(fracPerf(:,i) > th, 1);
% end

subplot(2,2,2)
median(xat(id))
histogram(xat(id), 'binEdges', 0:.1:1, 'Normalization', 'probability', 'FaceColor', .5*[1 1 1])
% xlabel(sprintf('Time until %d%% performance', round(100*th)))
%%

subplot(2,2,3)

x   = [Stmp(:).decodingAccuracy];
xci = [Stmp(:).decodingAccuracyCI];
y   = [Sfull(:).decodingAccuracy];
yci = [Sfull(:).decodingAccuracyCI];
plot([x; x], [y-yci; y+yci], 'Color', .5*[1 1 1]); hold on
plot([x-xci; x+xci], [y; y], 'Color', .5*[1 1 1]);
plot(x,y, 'k.')
hold on
xd = [.6 1];
plot(xd, xd, 'k')
xlim(xd); ylim(xd)

[~, pval] = ttest(x,y);

signrank(x,y)

y-x
%%
X = cell2mat(arrayfun(@(x,y) x.stim.pulses(y.trialIdx,:), O(:), Sfull(:), 'UniformOutput', false));
% y = cell2mat(arrayfun(@(x) x.yhat==1, Sfull(:), 'UniformOutput', false));
y = cell2mat(arrayfun(@(x) x.yhat==1, Stmp(:), 'UniformOutput', false));


[b, ~, stats] = glmfit(X,y, 'binomial');

ptime = 1:7;
subplot(2,2,4)
plot(ptime, b, 'k-')
hold on
plot(ptime, b, 'k.', 'MarkerSize', 10)
plot(ptime, b+stats.se, 'k--')
plot(ptime, b-stats.se, 'k--')
xlabel('Pulse #')
ylabel('Psychophysical Weight')
xlim([.5 7.5])
set(ax, 'Box', 'off', 'XTick', ptime, 'YTick', 0:.2:.8)

%%


% y = cell2mat(arrayfun(@(x) x.yhat==1, Sb(:), 'UniformOutput', false));
% y = cell2mat(arrayfun(@(x) x.yhat==1, Stmp(:), 'UniformOutput', false));


[b, ~, stats] = glmfit(X,y, 'binomial');

ptime = 1:7;
ax = subplot(2,2,4);
plot(ptime, b, 'r-')
hold on
plot(ptime, b, 'r.', 'MarkerSize', 10)
plot(ptime, b+stats.se, 'r--')
plot(ptime, b-stats.se, 'r--')
xlabel('Pulse #')
ylabel('Psychophysical Weight')
xlim([.5 7.5])
set(ax, 'Box', 'off', 'XTick', ptime, 'YTick', 0:.2:.8)


% %% Test whether quadratic terms help
% clear S1 S0
% for kEx = 1:numel(D)
%     
% [R, C] = D(kEx).getPredictors('quadratic', true, 'zscore', true);
% 
% S1(kEx) = decodeFun(R, C, decodeOpts{:});
% 
% [R, C] = D(kEx).getPredictors('quadratic', false, 'zscore', true);
% 
% S0(kEx) = decodeFun(R, C, decodeOpts{:});
% 
% fprintf('%s: %02.2f, %02.2f\n', D(kEx).exname, S0(kEx).decodingAccuracy, S1(kEx).decodingAccuracy)
% end

%%
w  = [];
pd = [];
nBins = min(arrayfun(@(x) numel(x.binIdx), Sfull));
for kEx = 1:numel(Stmp)
    neuronix = Stmp(kEx).wts~=0;
    pd_ = D(kEx).prefDir(neuronix);
    pd = [pd pd_(:)'];
    w_ = reshape(Sfull(kEx).wts, [], sum(neuronix));
    w  = [w w_(1:nBins,:)];
end

%%

w(:,isnan(pd)) = [];
pd(isnan(pd)) = [];
[~, id] = sort(pd);

[xx,yy] = meshgrid(linspace(-5,5, 10));
gkern = exp( - (xx.^2 + yy.^2));
I = conv2(w(:,id), gkern, 'same');
imagesc(I); colormap gray

%%
% 
% imagesc(xat(1:nBins), pd(id), I')

bs = 30;
binEdge = -180:bs:180;
[cnt, id] = histc(pd, binEdge);

bs = 30;
binEdge = 0:bs:180;
[~, id] = histc(pd, binEdge);

clf
ax = subplot(1,1,1);
n = numel(unique(id));
I = zeros(n, nBins);
cmap = pdsa.cbrewer('jake', 'rdbu', n);

tt = O(1).spikeBinCenters(Sfull(1).binIdx);

for i = 1:n
   iix = id==i;
   I(i,:) = mean(w(:,iix),2)';
   
   plot(tt, I(i,:), 'Color', cmap(i,:)); hold on
end




%%


binCenter = binEdge(1:end-1)+bs/2;
imagesc(xat(1:nBins), binCenter, I)
axis xy
set(ax, 'YTick', -180:90:180)
colormap jet
xlabel('Time')
ylabel('Preferred Direction')


%%
wtmp = bsxfun(@rdivide, w, sum(w));
wtmp = bsxfun(@times, w, sign(sum(w)));
clf
xat = D(1).spikeBinCenters(Sfull(1).binIdx);
plot(xat(1:nBins), wtmp); hold on
%%
clf
m = mean(wtmp,2);
s = std(wtmp, [],2)/sqrt(size(wtmp,2));
stairs(xat(1:nBins)-mean(diff(xat))/2, m, '-k'); hold on
stairs(xat(1:nBins)-mean(diff(xat))/2, m-s, '--');
stairs(xat(1:nBins)-mean(diff(xat))/2, m+s, '--');


% figt = figure_NmfPmfSummary(O, N, Sfull,1);

ppc = arrayfun(@(x) x.popnmf.pc, N);
mpc = arrayfun(@(x) x.pmf.pc, N);

ppc = arrayfun(@(x) x.popnmf.alpha, N);
mpc = arrayfun(@(x) x.pmf.alpha, N);

figure(1); clf
plot(mpc, ppc, '.'); hold on; plot(xlim, xlim,'k')

signrank(mpc, ppc)

%%
fig4(1) = figure_NmfPmfSummary(D, N, Sfull, 1);

pdsa.fixfigure(fig4(1), 8, [7.5 3])
saveas(fig4(1), fullfile(figDir, 'fig04_nmf-full-ver01.pdf'))

fig4(1) = figure_NmfPmfSummary(D, N, Sfull, 0);

pdsa.fixfigure(fig4(1), 8, [7.5 2.5])
saveas(fig4(1), fullfile(figDir, 'fig04_nmf-full-ver02.pdf'))

fig4(1) = figure_NmfPmfSummary(D, N, Sfull, 2);

pdsa.fixfigure(fig4(1), 8, [5.5 5.5])
saveas(fig4(1), fullfile(figDir, 'fig04_nmf-full-ver03.pdf'))

%%

fig5 = figure_CompareDecoders(D, Sfull, Sb);


%%
exampleSessions={'p20140304', 'p20140305', 'n20150316c', 'n20150324a'};
nExamples = numel(exampleSessions);
sx = ceil(sqrt(nExamples));
sy = round(sqrt(nExamples));

fig3 = figure(128); clf
clr = 'r';
for k = 1:nExamples
    
    subplot(sy,sx,k)
    
    kEx = find(strcmp(exampleSessions{k}, experiments));
    D(kEx).plotPMF(N(kEx).popnmf, 'r')
    D(kEx).plotPMF(N(kEx).pmf, 'k')
    [~,id]=min(N(kEx).nmf.alpha);
    plot(N(kEx).nmf.x, N(kEx).nmf.y(:,id), '.', 'Color', .5*[1 1 1])
    plot(N(kEx).nmf.xx(:,id), N(kEx).nmf.yy(:,id), 'Color', .5*[1 1 1])
    ylim([.4 1])
    xlabel('Motion Strength')
    ylabel('Proportion Correct')
    xlim([0 1.65])
    
    if kEx==1
        set(gca, 'Xtick', 0:.5:1.5)
        ylim([.5 1])
        text(.6,.65, 'Monkey')
        text(.6,.6, 'Population')
        text(.6,.55, 'Best Neuron')
        plot(.5+[0 .05], .65*[1 1], 'k')
        plot(.5+[0 .05], .6*[1 1], 'Color', clr)
        plot(.5+[0 .05], .55*[1 1], 'Color', .5*[1 1 1])
    end
    
    title('')
end

pdsa.fixfigure(fig3, 8, [5.5 5.5])
saveas(fig3, fullfile(figDir, 'fig03_examples-full-ver01.pdf'))