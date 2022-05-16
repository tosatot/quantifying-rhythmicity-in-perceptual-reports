function [tVal_Ob, tVal_Prm] = randEff(compSpectra)

% The function randEff computes a random effect statistical analysis

% Input Arguments:
%     compSpectra = three dimensional matrix containing complex values. The 3 dimensions
%     represent respectively iterations, frequencies, and participants.

% Output Arguments:
%     tVal_Ob     = t-values obtained by comparing the observed spectra to
%     the bias estimate spectra
%     tVal_Prm    = t-values obtained by comparing 2 groups to whom 
%     observed and bias estimate spectra are randomly assigned 
 
    
    totIt  = size(compSpectra,1);
    totFOI = size(compSpectra,2);
    totPr  = size(compSpectra,3);

    obs_spect     = squeeze(compSpectra(1,:,:));
    biasEst_spect = squeeze(mean(compSpectra(2:end,:,:),1));

    permMtrx      =randi(2,[totIt, totPr]);
    
    perm_obs     = nan(totIt,totFOI,totPr);
    perm_biasEst = nan(totIt,totFOI,totPr);

    for i =1:totIt
        perm_obs(i,:,:)     = cat(2,obs_spect(:,permMtrx(i,:)==1),biasEst_spect(:,permMtrx(i,:)==2));
        perm_biasEst(i,:,:) = cat(2,obs_spect(:,permMtrx(i,:)==2),biasEst_spect(:,permMtrx(i,:)==1));
    end

    tVal_Ob  =         abs(pairedTTest(obs_spect', biasEst_spect'));
    tVal_Prm = permute(abs(pairedTTest(permute(perm_obs,[3,2,1]),permute(perm_biasEst,[3,2,1]))),[3,2,1]);

end