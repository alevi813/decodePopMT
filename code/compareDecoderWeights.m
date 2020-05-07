choPath = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/residuals/data';
dirPath = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/residuals/data';

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

