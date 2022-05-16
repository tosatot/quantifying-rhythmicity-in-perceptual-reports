function out=sinFit(cfg,data,it)

% The function SinFit fits to the mean accuray time course a sinusoidal function 
% and a dampened harmonic oscillation. 

% Input Arguments:
%     cfg  = structure containing the parameters for the computation
%     data = structure containing the data 
%     it   = iteration number. when it > 1,  POI indixes are permuted to calculate statistics

% Output Arguments:
%     out = structure containing the cofficients of the fitted models, the
%     r-squared, and the RSS


    POI  = data.POI;    % probe onset interval
    BRV  = data.BRV;    % behavioural response value
    valid= data.valid;  % valid trials
    totPr= data.totPr;  % total number of participants
    
    TOI  = cfg.TOI;    % time of interest
    sigma= cfg.sigma;  % sigma for gaussian convolution
    dtOrd= cfg.dtOrd;  % DT order (reccomended to 1)

    
    atc_dt=nan(length(TOI),length(totPr));
    
    
    %______________calculate mean Accuracy Time Course ______________
    for pr=1:totPr

        POIsp = POI(valid(:,pr),pr);  %  POI of valid trials of a single participant 
        BRVsp = BRV(valid(:,pr),pr);  %  BRV of valid trials of a single participant 
      
        %SHUFFLE
        if it~=1   
            POIsp=POIsp(randperm(length(POIsp)));
        end

        %CONVOLUTION WITH A GAUSSIAN WINDOW
        W = exp(-((POIsp-TOI).^2)/(2*sigma^2));
        H = W.* BRVsp(:);
        Ht = nansum(H,1);
        Wt = nansum(W,1);
        ATC = Ht./Wt;

        %DETREND                                    % another option: accDt(:,pr) = detrend(acc(:,pr))
        [pCoef, s1,s2] = polyfit(TOI,ATC,dtOrd); % polynomial fit 
        trend = polyval(pCoef,TOI,s1,s2);
        atc_dt(:,pr) = ATC-trend;
        
    end
    
    mATC= mean(atc_dt,2);

    
    %__________________________fit SINE ______________________________
    [F0,gof]=fit(TOI(:),mATC(:),'sin1', 'Lower',[0,2*pi,-pi],'Upper',[2,2*pi*60,+pi]);
    

    %OUTPUTS #1
    out.acc    = mATC;        
    out.coeff  = coeffvalues(F0);
    out.rss    = gof.sse;   % RSS (Residual Sum of Squares) aka SSE (Sum of Squared Errors)
    out.rsquare= gof.rsquare;

        

    %___________fit Dampened Harmonic Oscillation _______________
    y  = mATC;
    x  = TOI(:);


    zc = x((diff(sign(y))~=0));         % Approximate Zero-Crossing Points
    yRng = max(y)- min(y);                % Estimate range of y
    yFrq = 1./(mean(diff(zc))*2);          % Estimate freq. (period multiplied by 2 to correct for the Zero-Crossings due to noise)
    yMn  = mean(y);                      % Estimate mean
    
    model= @(b,x)  b(1) .* exp(b(2).*x) .* (sin(2*pi*x.*b(3) + b(4))) + b(5); 
    [results, betas] = evalc('lsqcurvefit(model, [yRng; 0;  yFrq;  0;  yMn], x, y, [0.001; -20;  1;  -pi;  -0.5], [1; 20;  40;  +pi;  0.5])');
    yfit=model(betas,x);
    
    TSS = sum((y-mean(y)).^2);          % TSS (Total Sum of Squares) aka SST (Sum of Squares Total)
    RSS = sum((y(:)-yfit(:)).^2);       % RSS (Residual Sum of Squares) aka SSE (Sum of Squared Errors)
    rsquare = 1-RSS/TSS ;  
    
  
    %OUTPUTS #2
    out.coeffDMP  = betas;
    out.rsquareDMP= rsquare;
    out.rssDMP    = RSS;
   
end