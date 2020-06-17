 choPath = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/residuals/data';
 dirPath = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/residuals/data';

%%
condition = 'good';
experiments = getExperimentsAnd({condition});

for iExp = 1:length(experiments)
    
   dirSession = load([dirPath filesep experiments{iExp} '.mat']);
   choSession = load([choPath filesep experiments{iExp} '.mat']);
   
   nNeurons(iExp) = dirSession.nNeurons;
   if dirSession.nNeurons >1
       tmpr = corrcoef(dirSession.wAll, choSession.wAll);
       rval(iExp) = tmpr(1,2);
   else
       rval(iExp) = NaN;
   end
end

%%

dS = getDecodingData('all', 'nonstimulus', 'stimulus');
%%

plot(dS.flat.dir.w,  dS.flat.cho.slooW, 'o', 'color', [0 0 .6])
flat_r = corrcoef(dS.flat.dir.w,  dS.flat.cho.slooW); hold on

plot(dS.late.dir.w,  dS.late.cho.slooW, 'o', 'color', [.6 .6 0])
late_r = corrcoef(dS.late.dir.w,  dS.late.cho.slooW);

plot(dS.early.dir.w,  dS.early.cho.slooW, 'o', 'color', [.6 0 0])
early_r = corrcoef(dS.early.dir.w,  dS.early.cho.slooW);

axis square

%%
subplot(121)
plot(stim.flat.dir.w,  nonstim.flat.dir.w, 'o', 'color', [0 0 .6])
[flat.dir.r, flat.dir.p] = corrcoef(stim.flat.dir.w,  nonstim.flat.dir.w); hold on

plot(stim.late.dir.w,  nonstim.late.dir.w, 'o', 'color', [.6 .6 0])
[late.dir.r, late.dir.p] = corrcoef(stim.late.dir.w,  nonstim.late.dir.w);

plot(stim.early.dir.w,  nonstim.early.dir.w, 'o', 'color', [.6 0 0])
[early.dir.r, early.dir.p] = corrcoef(stim.early.dir.w,  nonstim.early.dir.w);

axis square

subplot(122)
plot(stim.flat.cho.slooW,  nonstim.flat.cho.slooW, 'o', 'color', [0 0 .6])
[flat.cho.r, flat.cho.p]    = corrcoef(stim.flat.cho.slooW,  nonstim.flat.cho.slooW); hold on

plot(stim.late.cho.slooW,  nonstim.late.cho.slooW, 'o', 'color', [.6 .6 0])
[late.cho.r, late.cho.p]   = corrcoef(stim.late.cho.slooW,  nonstim.late.cho.slooW);

plot(stim.early.cho.slooW,  nonstim.early.cho.slooW, 'o', 'color', [.6 0 0])
[early.cho.r, early.cho.p] = corrcoef(stim.early.cho.slooW,  nonstim.early.cho.slooW);

axis square


%%
dw1 = nonstim.flat.dir.w( nonstim.flat.dir.w  > -3e-3 & nonstim.flat.dir.w  < 3e-3 );
dw2 = stim.flat.dir.w( nonstim.flat.dir.w  > -3e-3 & nonstim.flat.dir.w  < 3e-3 );
[flat.dir.r, flat.dir.p] = corrcoef(dw1,  dw2);

cw1 = nonstim.late.cho.slooW( nonstim.late.cho.slooW  > -0.01 & nonstim.late.cho.slooW  < 0.01 );
cw2 = stim.late.cho.slooW( nonstim.late.cho.slooW  > -0.01 & nonstim.late.cho.slooW  < 0.01 );
[late.cho.r, late.cho.p] = corrcoef(cw1,  cw2); 

% flat dir p = 0.1195
% flat cho p = 