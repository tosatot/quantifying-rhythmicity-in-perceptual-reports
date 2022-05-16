function [data, vr] = genData(cfg, sd)

% The function 'genData' generates a set of simulated behavioural data with
% underlying rhythmicity

% Input Arguments:
%     cfg  = structure containing the parameters for the computation
%     sd   = scalar, which is the seed of the random number generator.

% Output Arguments:
%     data = contains the data generated in the simulation, to be used in the following analyses.
%     vr   = contains the distributions of variables used to generate the data


    rng(sd)  %seed random number generator
    
    % ____VARIABLES_____
    
    totPr      = cfg.totPr;
    totTr      = cfg.totTr;
    plvPr      = cfg.plvPr;
    plvTr      = cfg.plvTr;
    meanModDep = cfg.modDepth;
    
    %FREQUENCY
    %distribution of mean modulation freq. across participants
    meanFreqD  = 10;  %Hz  
    sigmaFreqD = 0.27;   
    %distribution of 'Sigmas' across participants
    meanSigmaD = 0.27; 
    sigmaSigmaD= 0.05;

    %PHASE
    %distribution of mean modulation phases across participants
    meanPhaseD = pi./2;
    kappaPhaseD= circ_kappa(plvPr);   %from PLV to K
    %distribution of 'Kappa' across participants
    meanKappaD = circ_kappa(plvTr);   %from PLV to K
    sigmaKappaD= 0.05;   
    

    %Filter
    bpWidth=0.2; %band pass width (+/-)
    order=3;     %order of the filter
    
    % Arousal 'weights'
    arousal=repmat(linspace(1,0,totTr)',1,totPr);
    
    %Initialize remaining varibles 
    modFreq  = nan(totTr,totPr);
    modPhase = nan(totTr,totPr);   
    POI      = nan(totTr,totPr);   %Probe Onset Interval
    BRV      = nan(totTr,totPr);   %Behavioural Response Values
    valid    = nan(totTr,totPr);
    prob     = nan(totTr,totPr,1000);
    
    
    %_____MAKE DISTRIBUTIONS____
    
    % Mean modulation Freq. across participant
    fd1 = makedist('Normal','mu',meanFreqD,'sigma',sigmaFreqD);
    meanFreq=random(fd1,[1,totPr]);    
    
    % Sigma (of mean mod. freq. ) across participant
    sd1 = makedist('Normal','mu',meanSigmaD,'sigma',sigmaSigmaD);
    sigmaFreq=random(sd1,[1,totPr]);
    sigmaFreq(sigmaFreq<0.01)=0.01;    %truncate distribution
    
    
    % Mean modulation Phase across participants
    if ~isinf(kappaPhaseD)   % 0 < plv < 1 
        meanPhase= circ_vmrnd(meanPhaseD,kappaPhaseD,totPr);   % Randomly sample a von Mises distribution, with preferred direction "meanPhasePr" and concentration parameter "kappaPhasesD".
    else    % plv = 1 
        meanPhase= ones(1,totPr).*meanPhaseD;                  % Participants are perfectly phase aligned
    end
    
    
    % Kappa (for mod. phases dist.) across participants
    kd1 = makedist('Normal','mu',meanKappaD,'sigma',sigmaKappaD);
    kappaTr=random(kd1,[1,totPr]); 
    
    
    % Mean Mod Depth 
    mdd=makedist('Normal','mu',meanModDep,'sigma',0.04);  
    modDepth=random(mdd,[1,totPr]); 
    if sum(modDepth > 0.49)>0
        modDepth(modDepth > 0.49)=0.49;
    end
    
    % Probability for valid trials
    pv1 = makedist('Normal','mu',0.80,'sigma',0.07);   %only 20% of trials is valid, so that the sampling (of the 1000 possible P.O.I.) is irregular 
    probValid=random(pv1,[1,totPr]);
    probValid(probValid>1)=1;
    probValid(probValid<0.05)=0.05;
    

    %_____ GENERATE DATA ______
    
    for pr=1:totPr    %_____Participants______
        
        POI(:,pr)     =  round(rand(totTr,1),3);     
        POI(POI<0.001)=  0.001;
        POIidx        =  POI(:,pr).*1000;            %round(linspace (1,1000, totTr)); 
 
        %frequency for each trial (Gaussian distribution)   
        fd2 = makedist('Normal','mu',meanFreq(pr),'sigma',sigmaFreq(pr));
        modFreq(:,pr)=random(fd2,[totTr,1]);
        
        %phase for each trial (Von Mises distribution)    
        modPhase(:,pr)= circ_vmrnd(meanPhase(pr),kappaTr(pr),totTr);

        % choose randomly which trials are valid 
        valid(:,pr)=datasample([0 1],totTr, 'Weights', [1-probValid(pr) probValid(pr)]);  
        
        
        for tr=1:totTr   %______Trials________
            
            %random noise
            tc=rand(32000,1);
            tc=tc-mean(tc);

            %filter 
            fs=1000; 
            fcutlow=modFreq(tr,pr)-bpWidth;   %low cut frequency in Hz
            fcuthigh=modFreq(tr,pr)+bpWidth;   %high cut frequency in Hz
            [b,a]=butter(order,[fcutlow,fcuthigh]/(fs/2));          
            f_tc_l=filter(b,a,tc);  
            
            %exclude "burn in" period
            f_tc=f_tc_l(30001:32000);  %filtered time course
            
            %Hilbert
            H1=hilbert(f_tc(1:900));
            [~, idx]=( min( abs( angle(H1(301:600)) - modPhase(tr,pr) ) ) ); %select starting phase


            filtSig=f_tc(idx+300:idx+300+999);  %select 1000 samples of filtered signal

            prob(tr,pr,:)=filtSig./max(abs(filtSig)).*modDepth(pr).*arousal(tr,pr)+0.5;   % the max amplitude of the signal will be +/- ( modulation depth of the sinusoidal + modulation of arousal)       
            
            prPOT=prob(tr,pr,POIidx(tr));  %find the probability for the BRVs of a given POT  
            BRV(tr,pr)=datasample ([-1 1],1, 'Weights', [1-prPOT prPOT]);   %calculate the BRV given the probability
        end
    end
    
    
    %_____ OUTPUT ________
    
    vr.prob     =prob;
    vr.meanFreq =meanFreq;
    vr.sigmaTr  =sigmaFreq;
    vr.modFreq  =modFreq;
    vr.meanPhase=meanPhase;
    vr.kappaTr  =kappaTr;
    vr.modPhase =modPhase;
    vr.modDepth =modDepth;
    
    data.POI     = POI;
    data.BRV     = BRV;
    data.valid   = logical(valid);
    data.arousal = arousal;
    data.totPr   = totPr;
end