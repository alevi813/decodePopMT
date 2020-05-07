%%
cptm = mean(dS.late.cho.cpt);

ptimes = 60:15:165;

p(1) = mean(cptm(ptimes(1):ptimes(2)));
p(2) = mean(cptm(ptimes(2):ptimes(3)));
p(3) = mean(cptm(ptimes(3):ptimes(4)));
p(4) = mean(cptm(ptimes(4):ptimes(5)));
p(5) = mean(cptm(ptimes(5):ptimes(6)));
p(6) = mean(cptm(ptimes(6):ptimes(7)));
p(7) = mean(cptm(ptimes(7):ptimes(8)));

%plot(1:7, p, '-o', 'color', [.6 0 0])

%% raw data
ptimes = 111:15:231;

p(1) = mean(cpm(ptimes(1):ptimes(2)));
p(2) = mean(cpm(ptimes(2):ptimes(3)));
p(3) = mean(cpm(ptimes(3):ptimes(4)));
p(4) = mean(cpm(ptimes(4):ptimes(5)));
p(5) = mean(cpm(ptimes(5):ptimes(6)));
p(6) = mean(cpm(ptimes(6):ptimes(7)));
p(7) = mean(cpm(ptimes(7):ptimes(8)));

plot(1:7, p, '-o', 'color', [.6 0 0])

% cp
% early: r= -0.816, p = 0.025
% flat: r= -0.037, p = 0.937
% late: r= -0.394, p = 0.381
%
% pta
% early: rval = 0.617, pval = 0.099
% flat:  rval = 0.839, pval = 0.019
% late:  rval = -0.734, pval = 0.060
%%
% early slope =  0.01437 +/- 0.0030 ci
% flat slope  =  0.01359 +/- 0.0022
% late slope  = 