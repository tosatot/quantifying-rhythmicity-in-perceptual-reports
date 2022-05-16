function out=stLSS(cfg,data,it)

% The function singTr applay the LSS and the WLSS methods to single trials. 

% Input Arguments:
%     cfg  = structure containing the parameters for the computation
%     data = structure containing the data 
%     it   = iteration number. when it > 1,  POI indixes are permuted to calculate statistics

% Output Arguments:
%     out = structure containing the complex spectra


    FOI    = cfg.FOI;    % frequency of interest   
    TOI    = cfg.TOI;    % time of interest
    sigma  = cfg.sigma;  % sigma for gaussian convolution
    dtOrd  = cfg.dtOrd;  % DT order (reccomended to 1)
    totPr  = cfg.totPr;  % total number of participants 
    
    POI    = data.POI;    % probe onset interval
    BRV    = data.BRV;    % behavioural response value
    valid  = data.valid;  % valid trials
    arousal= data.arousal;% arousal, or any mesured variable which may influence the strengh of the rhythm 
  
    stLSS   = nan(length(FOI),totPr);
    stWLSS  = nan(length(FOI),totPr);

    
    for pr =1:totPr
        
        POIsp = POI(valid(:,pr),pr);  %  POI of valid trials of a single participant 
        BRVsp = BRV(valid(:,pr),pr);  %  BRV of valid trials of a single participant 
        
        ar  = arousal(valid(:,pr),pr)./nanmean(arousal(valid(:,pr),pr));  %normalize arousal values, so that they integrate to 1

        %SHUFFLE
        if it~=1   
            POIsp=POIsp(randperm(length(POIsp)));
        end
        
        %DETREND
        %CONVOLUTION WITH GAUSSIAN ON THE TIME DOMAIN
        W = exp(-((POIsp-TOI).^2)/(2*sigma^2));
        H = W.* BRVsp(:);
        Ht = nansum(H,1);
        Wt = nansum(W,1);
        atc= Ht./Wt;            % accuracy time course
        %SUBTRACTION OF THE TREND TO SINGLE TRIALS
        [pCoef, s1,s2] = polyfit(TOI,atc,dtOrd); %polynomial fit 
        trend = polyval(pCoef,TOI,s1,s2);
        BRVsp=BRVsp-trend(POIsp.*1000)';
        
        %TAPERING OF SINGLE TRIALS
        taper=hann(1000);
        taper=taper./sum(taper);
        BRVsp=BRVsp.*(taper(POIsp.*1000));
        
        %DEMEAN again
        BRVsp=BRVsp-mean(BRVsp);
        
        
        
        %____ stLSS ____
        for fr = 1:length(FOI)
            POP = POIsp*2*pi*FOI(fr);  %Probe Onset Phase
            sinePOP = sin(POP);
            cosPOP  = cos(POP);    
            
            X  = [ones(length(sinePOP),1), cosPOP, sinePOP]; 
            Y  = BRVsp;
            betas = X\Y;
            
            stLSS(fr,pr) = complex(betas(2),betas(3));  
        end

        
        %____ stWLSS ____
        for fr = 1:length(FOI)
            POP = POIsp*2*pi*FOI(fr);
            sinePOP = sin(POP);
            cosPOP  = cos(POP);  
            
            X  = [ones(length(POIsp),1), cosPOP.*ar, sinePOP.*ar]; 
            Y  = BRVsp;
            betas = X\Y;
            
            stWLSS(fr,pr) = complex(betas(2),betas(3));
        end
        
    end


    out.stLSS =stLSS;
    out.stWLSS=stWLSS;
end

