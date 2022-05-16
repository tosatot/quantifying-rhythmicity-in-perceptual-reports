function out=atcDFT(cfg,data,it)

% The function atcDFT applay the DFT to the mean accuracy time course.

% Input Arguments:
%     cfg  = structure containing the parameters for the computation
%     data = structure containing the data 
%     it   = iteration number. when it > 1,  POI indixes are permuted to calculate statistics

% Output Arguments:
%     out = structure containing the complex spectra


    FOI  = cfg.FOI;    % frequency of interest   
    TOI  = cfg.TOI;    % time of interest
    sigma= cfg.sigma;  % sigma for gaussian convolution
    dtOrd= cfg.dtOrd;  % DT order (reccomended to 1)
    totPr= cfg.totPr;  % total number of participants
    
    POI  = data.POI;    % probe onset interval
    BRV  = data.BRV;    % behavioural response value
    valid= data.valid;  % valid trials

    atcDFT = nan(length(FOI),totPr);

    
    for pr= 1:totPr

        %LOAD 
        POIsp = POI(valid(:,pr),pr);  %  POI of valid trials of a single participant 
        BRVsp = BRV(valid(:,pr),pr);  %  BRV of valid trials of a single participant 

        %SHUFFLE
        if it~=1   
            POIsp=POIsp(randperm(length(POIsp)));  %Suffle the POTs to calculate permutation statistics (for all the iterations different from 1, the first iteration is the "obserbed data") 
        end
        
        %CONVOLUTION WITH A GAUSSIAN WINDOW
        W = exp(-((POIsp-TOI).^2)/(2*sigma^2));
        H = W.* BRVsp(:);
        Ht = nansum(H,1);
        Wt = nansum(W,1);
        atc = Ht./Wt;
       
        %DETREND                                             
        [pCoef, s1,s2] = polyfit(TOI,atc,dtOrd);       
        trend = polyval(pCoef,TOI,s1,s2);
        atc_dt = atc-trend;
         

        %DFT OF MEAN ACCURACY TIME COURSE
        fs = 1000;
        Nt = length(TOI);

        padTo = 4000;
        nHz = floor(padTo/2)+1;
        hz = linspace(0,0.5,nHz)*fs;

        taper = hann(Nt);
        taper = taper./sum(taper);

        xfft = fft(atc_dt(:).*taper,padTo).*2;
        xfft = xfft(1:nHz);
        
        atcDFT(:,pr)  = xfft(ismembertol(hz,FOI,0.0001),:);
        
 
    end
    
    out.atcDFT = atcDFT;

end

