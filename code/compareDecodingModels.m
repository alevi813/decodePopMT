comp = getComp;

if strcmp(comp, 'laptop')
    choDir = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/choice/data';
    dirDir = '/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/data';
    
    exNames = dir('/Users/aaronlevi/Dropbox/twagAnalysis4.1/decoding/direction/data/');
else
    choDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/choice/data';
    dirDir = '/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/direction/data';
    
    exNames = dir('/Users/Aaron/Dropbox/twagAnalysis4.1/decoding/direction/data/');
end

exNames = exNames(4:end);

wAll = nan(length(exNames), 2);

dvDistance.across     = nan(length(exNames), 1);
dvThetaDegrees.across = nan(length(exNames), 1);
wDistance.across      = nan(length(exNames), 1);
wThetaDegrees.across  = nan(length(exNames), 1);

dvDistance.within     = nan(length(exNames), 1);
dvThetaDegrees.within = nan(length(exNames), 1);
wDistance.within      = nan(length(exNames), 1);
wThetaDegrees.within  = nan(length(exNames), 1);

nNeurons       = nan(length(exNames), 1);

for kExp = 1:length(exNames)
%for kExp = 2
    
    % get the ensemble activity for choice and direction models for a session
    dirSession = load([dirDir filesep exNames(kExp).name]);
    choSession = load([choDir filesep exNames(kExp).name]);
    
%% new regression
%  response (1st column dir spikes, 2nd column cho spikes) and choice
    popR = [sum(dirSession.dvAll, 2) , sum(choSession.dvAll, 2)]; 
    w2  = glmfit(popR,  dirSession.cho); 

    allW(kExp, :) = w2(2:end);


%%  dvs   
    % ACROSS
    % take the mean, shorten it to get rid of peri-saccadic activity
    mdvDirAll = mean(dirSession.dvAll); mdvDirAll = mdvDirAll(1:120);
    mdvChoAll = mean(choSession.dvAll); mdvChoAll = mdvChoAll(1:120);
    
    % get euclidean distance of the two mean dvs
    dvDistance.across(kExp) = norm(mdvDirAll - mdvChoAll);
    
    % get norms
    normDir = norm(mdvDirAll);
    normCho = norm(mdvChoAll);
    
    % get cosine similarity
    % dot product of the two model dvs, normalized by euclidean norm
    % MEAN DVs
    dvCosTheta = dot(mdvDirAll, mdvChoAll) / (normDir*normCho);
    dvThetaDegrees.across(kExp) = acosd(dvCosTheta);
    
    % WITHIN
    % take the mean, shorten it to get rid of peri-saccadic activity
    mdvTrain = mean(dirSession.dvTrain); mdvTrain = mdvTrain(1:120);
    mdvTest  = mean(dirSession.dvTest); mdvTest = mdvTest(1:120);
    
    % get euclidean distance of the two mean dvs
    dvDistance.within(kExp) = norm(mdvTrain - mdvTest);
    
    % get norms
    normTrain = norm(mdvTrain);
    normTest  = norm(mdvTest);
    
    % get cosine similarity
    % dot product of the two model dvs, normalized by euclidean norm
    % MEAN DVs
    dvCosTheta = dot(mdvTrain, mdvTest) / (normTrain*normTest);
    dvThetaDegrees.within(kExp) = acosd(dvCosTheta);    
    
%%    
    % dv by trial
    dvDir = dirSession.dvAll(:,1:120);
    dvCho = choSession.dvAll(:,1:120);    
    
    dvDirTrain = dirSession.dvTrain(:,1:120);
    dvChoTrain = choSession.dvTrain(:,1:120);  
    
    dvDirTest = dirSession.dvTest(:,1:120);
    dvChoTest = choSession.dvTest(:,1:120);  

%% over time    
    % get cosine similarity
    % dot product of the two model dvs, normalized by euclidean norm
    thetaTimeAll   = nan(1,120);
    thetaTimeTrain = nan(1,120);
    thetaTimeTest  = nan(1,120);
    
    for iBin = 1:120
        % all data ACROSS dir/cho
        tempDot     = dot(dvDir(iBin,:), dvCho(iBin,:));
        tempNormDir = norm(dvDir(iBin,:));
        tempNormCho = norm(dvCho(iBin,:));
        
        thetaTime.across(iBin) = acosd( tempDot / (tempNormDir*tempNormCho) );
        
        % WITHIN dir model   
        tempDot       = dot(dvDirTrain(iBin,:), dvDirTest(iBin,:));
        tempNormTrain = norm(dvDirTrain(iBin,:));
        tempNormTest  = norm(dvDirTest(iBin,:));
        
        thetaTime.within(iBin) = acosd( tempDot / (tempNormDir*tempNormTest) );  
        
    end

%%  compare weights
    % get euclidean distance of the weights themselves
    % ACROSS models
    wDistance.across(kExp) = norm(dirSession.wAll - choSession.wAll);
    % WITHIN (direction) model
    wDistance.within(kExp) = norm(dirSession.wTrain - dirSession.wTest);
    
    % cosine similarity of weights
    % ACROSS models
    normDirW = norm(dirSession.wAll);
    normChoW = norm(choSession.wAll);
        
    wCosTheta = dot(dirSession.wAll, choSession.wAll) / (normDirW*normChoW);
    wThetaDegrees.across(kExp) = acosd(wCosTheta);    
    
    % WITHIN models
    normWtrain = norm(dirSession.wTrain);
    normWtest = norm(dirSession.wTest);
        
    wCosTheta = dot(dirSession.wTrain, dirSession.wTest) / (normWtrain*normWtest);
    wThetaDegrees.within(kExp) = acosd(wCosTheta);
    
    % how many neurons
    nNeurons(kExp) = length(dirSession.wAll);
end % experiment loop

%% plot some stuff (across)
% dvDistance.across     
% dvThetaDegrees.across
% wDistance.across     
% wThetaDegrees.across  


figure
plot(wDistance.across, 'o')
hold on
plot(dvDistance.across, 'o')

xlabel('session')
ylabel('euclidean distance')
legend('weights', 'ensemble activity')

figure
wThetaDegrees.across(imag(wThetaDegrees.across) ~=0) = 0;
plot(wThetaDegrees.across, 'o')

hold on
dvThetaDegrees.across(imag(dvThetaDegrees.across) ~=0) = 0;
plot(dvThetaDegrees.across, 'o')

ylabel('angle (degrees)')
xlabel('session')
legend('weights', 'ensemble activity')

yticks([0:30:180]);

%% plot some stuff (within)
figure
plot(wDistance.within, 'o')
hold on
plot(dvDistance.within, 'o')

xlabel('session')
ylabel('euclidean distance')
legend('weights', 'ensemble activity')

wThetaDegrees.within(imag(wThetaDegrees.within) ~=0) = 0;
dvThetaDegrees.within(imag(dvThetaDegrees.within) ~=0) = 0;

figure
plot(wThetaDegrees.within, 'o')
hold on
plot(dvThetaDegrees.within, 'o')

ylabel('angle (degrees)')
xlabel('session')
legend('weights', 'ensemble activity')

yticks([0:30:180]);
%% test
pWtheta  = ranksum(wThetaDegrees.within, wThetaDegrees.across);
pDVtheta = ranksum(dvThetaDegrees.within, dvThetaDegrees.across);

pWdistance  = ranksum(wDistance.within, wDistance.across);
pDVdistance = ranksum(dvDistance.within, dvDistance.across);
