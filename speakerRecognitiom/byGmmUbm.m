%
%
%  Speaker Recognition based on GMM UBM
%  This program is based on Matlab MSR toolbox
%  Created on: Feb 12, 2016
%  Author: Adam
%

nSpeakers = 20;  
nDims = 13;             % dimensionality of feature vectors  
nMixtures = 32;         % How many mixtures used to generate data  
nChannels = 10;         % Number of channels (sessions) per speaker  
nFrames = 100;         % Frames per speaker (10 seconds assuming 100 Hz)  
nWorkers = 1;           % Number of parfor workers, if available  
% Pick random centers for all the mixtures.  
mixtureVariance = .10;  
channelVariance = .05;  
mixtureCenters = randn(nDims, nMixtures, nSpeakers);  
channelCenters = randn(nDims, nMixtures, nSpeakers, nChannels)*.1;  
trainSpeakerData = cell(nSpeakers, nChannels);  
testSpeakerData = cell(nSpeakers, nChannels);  
speakerID = zeros(nSpeakers, nChannels);  
  
% Create the random data. Both training and testing data have the same  
% layout.  
disp('Create the random data');  
for s=1:nSpeakers  
    trainSpeechData = zeros(nDims, nMixtures);  
    testSpeechData = zeros(nDims, nMixtures);  
    for c=1:nChannels  
        for m=1:nMixtures  
            % Create data from mixture m for speaker s  
            frameIndices = m:nMixtures:nFrames;  
            nMixFrames = length(frameIndices);  
            trainSpeechData(:,frameIndices) = ...  
                randn(nDims, nMixFrames)*sqrt(mixtureVariance) + ...  
                repmat(mixtureCenters(:,m,s),1,nMixFrames) + ...  
                repmat(channelCenters(:,m,s,c),1,nMixFrames);  
            testSpeechData(:,frameIndices) = ...  
                randn(nDims, nMixFrames)*sqrt(mixtureVariance) + ...  
                repmat(mixtureCenters(:,m,s),1,nMixFrames) + ...  
                repmat(channelCenters(:,m,s,c),1,nMixFrames);  
        end  
        trainSpeakerData{s, c} = trainSpeechData;  
        testSpeakerData{s, c} = testSpeechData;  
        speakerID(s,c) = s;                 % Keep track of who this is  
    end  
end  

% Step1: Create the universal background model from all the training speaker data  
disp('Create the universal background model');  
nmix = nMixtures;           % In this case, we know the # of mixtures needed  
final_niter = 10;  
ds_factor = 1;  
ubm = gmm_em(trainSpeakerData(:), nmix, final_niter, ds_factor, nWorkers);  
% Step2: Now adapt the UBM to each speaker to create GMM speaker model.  
disp('Adapt the UBM to each speaker');  
map_tau = 10.0;  
config = 'mwv';  
gmm = cell(nSpeakers, 1);  
for s=1:nSpeakers  
    disp(['for the ',num2str(s),' speaker...']);  
    gmm{s} = mapAdapt(trainSpeakerData(s, :), ubm, map_tau, config);  
end  
% Step3: Now calculate the score for each model versus each speaker's data.  
% Generate a list that tests each model (first column) against all the  
% testSpeakerData.  
trials = zeros(nSpeakers*nChannels*nSpeakers, 2);  
answers = zeros(nSpeakers*nChannels*nSpeakers, 1);  
for ix = 1 : nSpeakers,  
    b = (ix-1)*nSpeakers*nChannels + 1;  
    e = b + nSpeakers*nChannels - 1;  
    trials(b:e, :)  = [ix% ones(nSpeakers*nChannels, 1), (1:nSpeakers*nChannels)'];  
    answers((ix-1)*nChannels+b : (ix-1)*nChannels+b+nChannels-1) = 1;  
end  
disp('Calculate the score for each model vs test speaker');  
gmmScores = score_gmm_trials(gmm, reshape(testSpeakerData', nSpeakers*nChannels,1), trials, ubm);  
% Step4: Now compute the EER and plot the DET curve and confusion matrix  
imagesc(reshape(gmmScores,nSpeakers*nChannels, nSpeakers))  
title('Speaker Verification Likelihood (GMM Model)');  
ylabel('Test # (Channel x Speaker)'); xlabel('Model #');  
colorbar; drawnow; axis xy  
figure  
disp('Compute the EER');  
[eer,auc] = compute_eer(gmmScores, answers, true);  