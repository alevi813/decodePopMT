% direction
%session = load('/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/direction/data/l20190806.mat');

% choice 
session = load('/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/choice/data/l20190806.mat');

spikeLeft  = session.dvAll(session.direction==-1,:);
spikeRight = session.dvAll(session.direction==1,:);

spikeLeft  = sum(spikeLeft(:, 50:170),2);
spikeRight = sum(spikeRight(:, 50:170),2);

spikesAll = sum(session.dvAll(:, 50:170), 2);

pulsecoh = sum(session.stim.pulses, 3);
sumcoh   = mean(pulsecoh, 2) ./19;

if isfield(session.stim, 'validTrials')
    pulsecoh = pulsecoh(session.stim.validTrials, :);
    sumcoh = sumcoh(session.stim.validTrials);
else
    pulsecoh = pulsecoh(session.stim.goodtrial, :);
    sumcoh = sumcoh(session.stim.goodtrial);    
end

neuralChoice = spikesAll>0;
%neuralChoice = spikesAll> -0.1355;

%%
figure
subplot(121)
pmf  = pmfTools(sumcoh, neuralChoice, 'nBins', 18);
fit(pmf);
plot(pmf);

subplot(122)
ppk = ppkTools(pulsecoh, neuralChoice, 'ridge', true);
plot(ppk, 'plotFit', false);

