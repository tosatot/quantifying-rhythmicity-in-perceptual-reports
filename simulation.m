%% add paths


% this script will guide you through the generation of simulated data, the
% analysis of those data, and the visualization.

% function from the "Circular Statistics Toolbox" are needed. see:
% P. Berens, CircStat: A Matlab Toolbox for Circular Statistics, Journal of Statistical Software, Volume 31, Issue 10, 2009

%% GENERATE DATA

cfg.totPr    = 30;    % total participannts
cfg.totTr    = 400;   % total trials per participant
cfg.plvTr    = 0.92;  % phase locking value between trials of the same participant
cfg.plvPr    = 0.92;  % phase locking value between participants
cfg.modDepth = 0.32;  % mean modulation depth

sd   = randi(1000,1); % seed of random number generator 

[data, vr] = genData(cfg,sd);

%% analyse simulated data

cfg.TOI   = 0.001: 0.001:1;  % Times Of Interest (sec.)
cfg.FOI   = 1:0.25:40;       % Frequences Of Interest (Hz)
cfg.sigma = 0.01;            % sigma of the gaussian for the convolution in the time domain (sec.)
cfg.dtOrd = 1; 
it=1;


%sine fitting
outSF=sinFit(cfg,data,it);

%mean accuracy time course methods
outATC=atcDFT(cfg,data,it);

%single trials methods
outST=stLSS(cfg,data,it);

%% calculate STATISTICS

totIt=1001; % total iterations. in the first iteration obseved data are analysed, in the other iterations the BRV vector is permuted. 

ZatcDFT = nan(totIt, length(cfg.FOI), data.totPr);
ZstLSS  = nan(totIt, length(cfg.FOI), data.totPr);
ZstWLSS = nan(totIt, length(cfg.FOI), data.totPr);
    
coeff       = nan(totIt, 3);
rsquare     = nan(totIt, 1);
coeffDMP    = nan(totIt, 5);
rsquareDMP  = nan(totIt, 1);
atc         = nan(totIt, 1000);

sd   = randi(1000,1); 
            
%___SPECTRAL DECOMPOSITION___
for it =1:totIt

    ot=sinFit(cfg,data,it);
    coeff(it,:)     = ot.coeff(1:3);
    rsquare(it,:)   = ot.rsquare;
    coeffDMP(it,:)  = ot.coeffDMP;
    rsquareDMP(it,:)= ot.rsquareDMP;
    atc(it,:)       = ot.acc;
    
    ot=atcDFT(cfg,data,it);
    ZatcDFT(it,:,:) = ot.atcDFT; 

    ot=stLSS(cfg,data,it);
    ZstLSS(it,:,:)  = ot.stLSS; 
    ZstWLSS(it,:,:) = ot.stWLSS; 
    
end

% For FIXED EFFECT
PatcDFT   = squeeze(abs(mean(ZatcDFT,3))).^2;     % power atcFFT 
PstLSS    = squeeze(abs(mean(ZstLSS,3))).^2;      % power stLSS  
PstWLSS   = squeeze(abs(mean(ZstWLSS,3))).^2;     % power stWLSS


% For RANDOM EFFECT
[tVatcDFT_Ob, tVatcDFT_Pr] = randEff(ZatcDFT);   % t-value atcFFT 
[tVstLSS_Ob,  tVstLSS_Pr]  = randEff(ZstLSS);    % t-value stLSS  
[tVstWLSS_Ob, tVstWLSS_Pr] = randEff(ZstWLSS);   % t-value stWLSS

%% PLOT Sine Fit

figure('position', [100 200 800 650],'Name',  ' sine fit ' )
set(gcf,'color','white')

subplot(2,2,1),
plot(cfg.TOI,atc(1,:),'lineWidth',1.4), hold on

plot(cfg.TOI, coeff(1,1)*sin(coeff(1,2)*cfg.TOI+coeff(1,3)),'lineWidth',1.4);
ylim([-0.4, 0.4])
xlabel('POI (ms)');
ylabel('Accuracy')
legend('Accuracy','Sine Fit')
legend boxoff
set(gca,'TickDir','out');
box off
text(-0.2,1.13,'A','units','normalized', 'FontSize', 16);
pos = get(gca, 'Position');
pos(4) = pos(4) -0.08;
pos(2) = pos(2) +0.02;
set(gca, 'Position', pos)

subplot(2,2,2),
histogram(rsquare(2:end),40,'facecolor',[0.2, 0.4, 0.2]), hold on, 
xline(rsquare(1),'r','lineWidth',1.5), 
xlim([0,rsquare(1).*1.2])
xlabel('R-squared')
ylabel('Nr. of Permutations')
set(gca,'TickDir','out');
box off
text(-0.2,1.13,'B','units','normalized', 'FontSize', 16);
pos = get(gca, 'Position');
pos(4) = pos(4)- 0.08;
pos(2) = pos(2) +0.02;
set(gca, 'Position', pos)

subplot(2,2,3),
plot(cfg.TOI,atc(1,:),'lineWidth',1.4), hold on
b=coeffDMP(1,:);
x=cfg.TOI;
plot(cfg.TOI, b(1) .* exp(b(2).*x) .* (sin(2*pi*x.*b(3) + b(4))) + b(5),'lineWidth',1.4);
ylim([-0.4, 0.4])
xlabel('POI (ms)');
ylabel('Accuracy')
legend('Accuracy','Damped Oscillation')
legend boxoff
set(gca,'TickDir','out');
box off
text(-0.2,1.13,'C','units','normalized', 'FontSize', 16);
pos = get(gca, 'Position');
pos(4) = pos(4) -0.08;
pos(2) = pos(2) +0.02;
set(gca, 'Position', pos)

subplot(2,2,4),
histogram(rsquareDMP(2:end),40,'facecolor',[0.2, 0.4, 0.2]), hold on, 
xline(rsquareDMP(1),'r','lineWidth',1.5), 
xlim([0, rsquareDMP(1).*1.2])
xlabel('R-squared')
ylabel('Nr. of Permutations')
set(gca,'TickDir','out');
box off
text(-0.2,1.13,'D','units','normalized', 'FontSize', 16);
pos = get(gca, 'Position');
pos(4) = pos(4)- 0.08;
pos(2) = pos(2) +0.02;
set(gca, 'Position', pos)

set(gcf,'Renderer', 'Painters'); 
% orient(gcf,'portrait')
% print(gcf,'Fig2- Sine Fitting','-dpdf')

%% PLOT SPECTRA (Multiple Comparison Methods & Fixed/Random Effect)  
figure('position',[20 0 801 1050],'Name',  ' freq. analysis and statistics ' ); %open figure 1
set(gcf,'color','white')

col1=[0.4940, 0.1840, 0.5560];
col2=[0.8500, 0.3250, 0.0980];
col3=[0.9290, 0.6940, 0.1250];

%______________Plot atcDFT______________________
subplot(6,3, 1), hold on; % F.E. bonferroni
average_pow=mean(PatcDFT(2:end,:),1);
plot(cfg.FOI,PatcDFT(1,:),'lineWidth',1.4);
plot(cfg.FOI,average_pow,'color',col3,'lineWidth',1.1);
plot(cfg.FOI,prctile(PatcDFT(2:end,:),95),'color',col2,'lineWidth',1.1); %alpha 0.7
plot(cfg.FOI,prctile(PatcDFT(2:end,:),99.9),'color',col1,'lineWidth',1.1);
set(gca,'TickDir','out');
set(gca,'YTick',[0, 0.01] )
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('Power');
xlabel('Freq. (Hz)');     
legend('signal','mean','95%','99.9%');
text(-0.3, 0.15,'Fix. Eff.','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.46,-0.55,'atc DFT','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.44,-1.2,'__________________________','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(0.3,1.2,'Bonferroni','units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.2,1.13,'A','units','normalized', 'FontSize', 16);

subplot(6,3,2), hold on;  % F.E. max based
plot(cfg.FOI,PatcDFT(1,:)./average_pow,'lineWidth',1.4);
yline(prctile(max(PatcDFT(2:end,:)./average_pow,[],2),95),'color',col1,'lineWidth',1.1,'alpha', 1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('Norm. Pow.');
xlabel('Freq. (Hz)');     
legend('signal','95%');
set(gca,'TickDir','out');
text(0.3,1.2,'Max-Based','units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.2,1.13,'B','units','normalized', 'FontSize', 16);

subplot(6,3,3), hold on,  %F.E. FDR
pval=(sum(PatcDFT(2:end,:)>PatcDFT(1,:),1)./(size(PatcDFT,1)-1));  %p-value  
fdr = calcFDR(pval);
[~,b]=find(fdr<0.1);
plot(cfg.FOI,PatcDFT(1,:),'lineWidth',1.4);
plot(cfg.FOI(b),ones(1,length(b)).*max(PatcDFT(1,:)).*1.1,'*','color',col1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('Power');
xlabel('Freq. (Hz)');  
legend('signal','sign. bin');
set(gca,'TickDir','out');
text(0.42,1.2,'FDR','units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.2,1.13,'C','units','normalized', 'FontSize', 16);

subplot(6,3,4), hold on  % R.E. bonferroni
plot(cfg.FOI, tVatcDFT_Ob,'lineWidth',1.4);
plot(cfg.FOI, mean(tVatcDFT_Pr,1),'color',col3,'lineWidth',1.1);
plot(cfg.FOI, prctile(tVatcDFT_Pr,95),'color',col2,'lineWidth',1.1);
plot(cfg.FOI, prctile(tVatcDFT_Pr,99.9),'color',col1,'lineWidth',1.1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('T-values');
xlabel('Freq. (Hz)');    
set(gca,'TickDir','out');
text(-0.3,0.15,'Rand. Eff.','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.2,1.13,'D','units','normalized', 'FontSize', 16);

subplot(6,3,5),hold on;  % R.E. max based
plot(cfg.FOI, tVatcDFT_Ob,'lineWidth',1.4);
yline(prctile(max(tVatcDFT_Pr,[],2),95),'color',col1,'lineWidth',1.1,'alpha', 1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('T-values');
xlabel('Freq. (Hz)');     
set(gca,'TickDir','out');
text(-0.2,1.13,'E','units','normalized', 'FontSize', 16);

subplot(6,3,6), hold on,  %F.E. FDR
pval=sum(tVatcDFT_Pr>tVatcDFT_Ob,1)./size(tVatcDFT_Pr,1);  %p-value  
fdr = calcFDR(pval);
[~,b]=find(fdr<0.1);
plot(cfg.FOI,tVatcDFT_Ob,'lineWidth',1.4);
plot(cfg.FOI(b),ones(1,length(b)).*max(tVatcDFT_Ob).*1.1,'*','color',col1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('T-values');
xlabel('Freq. (Hz)');  
set(gca,'TickDir','out');
text(-0.2,1.13,'F','units','normalized', 'FontSize', 16);


%________________Plot stLSS__________________
subplot(6,3,7), hold on;
plot(cfg.FOI, PstLSS(1,:),'lineWidth',1.4); % F.E. bonferroni
plot(cfg.FOI, mean(PstLSS(2:end,:),1),'color',col3,'lineWidth',1.1);
plot(cfg.FOI, prctile(PstLSS(2:end,:),95,1),'color',col2,'lineWidth',1.1);
plot(cfg.FOI, prctile(PstLSS(2:end,:),99.9,1),'color',col1,'lineWidth',1.1);
xlim([cfg.FOI(1) cfg.FOI(end)])
xlabel('Freq. (Hz)');
ylabel('Power');
set(gca,'TickDir','out');
text(-0.3, 0.15,'Fix. Eff.','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.46,-0.5,'st LSS','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.44,-1.2,'__________________________','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.2,1.13,'G','units','normalized', 'FontSize', 16);

subplot(6,3,8), hold on; % F.E. max space
plot(cfg.FOI, PstLSS(1,:),'lineWidth',1.4);
yline(prctile(max(PstLSS(2:end,:),[],2),95),'color',col1,'lineWidth',1.1,'alpha', 1);
xlim([cfg.FOI(1) cfg.FOI(end)])   
ylabel('Power');
xlabel('Freq. (Hz)');     
set(gca,'TickDir','out');
text(-0.2,1.13,'H','units','normalized', 'FontSize', 16);

subplot(6,3,9), hold on;%F.E. FDR
pval=sum(PstLSS(2:end,:)>PstLSS(1,:),1)./(size(PstLSS,1)-1);  %p-value  
fdr = calcFDR(pval);
[~,b]=find(fdr<0.1);
plot(cfg.FOI,PstLSS(1,:),'lineWidth',1.4);
plot(cfg.FOI(b),ones(1,length(b)).*max(PstLSS(1,:)).*1.1,'*','color',col1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('Power');
xlabel('Freq. (Hz)'); 
text(-0.2,1.13,'I','units','normalized', 'FontSize', 16);

subplot(6,3,10), hold on  % R.E. bonferroni
plot(cfg.FOI, tVstLSS_Ob,'lineWidth',1.4);
plot(cfg.FOI, mean(tVstLSS_Pr,1),'color',col3,'lineWidth',1.1);
plot(cfg.FOI, prctile(tVstLSS_Pr,95,1),'color',col2,'lineWidth',1.1);
plot(cfg.FOI, prctile(tVstLSS_Pr,99.9,1),'color',col1,'lineWidth',1.1);
xlim([cfg.FOI(1) cfg.FOI(end)])
ylabel('T-values')
xlabel('Freq. (Hz)')     
set(gca,'TickDir','out');
text(-0.3,0.15,'Rand. Eff.','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.2,1.13,'J','units','normalized', 'FontSize', 16);

subplot(6,3,11),hold on;  % R.E. max based
plot(cfg.FOI, tVstLSS_Ob,'lineWidth',1.4);
yline(prctile(max(tVstLSS_Pr,[],2),95),'color',col1,'lineWidth',1.1,'alpha', 1);
xlim([cfg.FOI(1) cfg.FOI(end)])
ylabel('T-values')
xlabel('Freq. (Hz)')     
set(gca,'TickDir','out');
text(-0.2,1.13,'K','units','normalized', 'FontSize', 16);

subplot(6,3,12), hold on,  %F.E. FDR
pval=sum(tVstLSS_Pr>tVstLSS_Ob,1)./size(tVstLSS_Pr,1);  %p-value  
fdr = calcFDR(pval);
[~,b]=find(fdr<0.1);
plot(cfg.FOI,tVstLSS_Ob,'lineWidth',1.4);
plot(cfg.FOI(b),ones(1,length(b)).*max(tVstLSS_Ob).*1.1,'*','color',col1);
xlim([cfg.FOI(1) cfg.FOI(end)])
ylabel('T-values')
xlabel('Freq. (Hz)') 
set(gca,'TickDir','out');
text(-0.2,1.13,'L','units','normalized', 'FontSize', 16);



%__________________Plot WLSS______________________
subplot(6,3,13), hold on; % F.E. bonferroni
plot(cfg.FOI, PstWLSS(1,:),'lineWidth',1.4);
plot(cfg.FOI,nanmean(PstWLSS(2:end,:),1),'color',col3,'lineWidth',1.1);
plot(cfg.FOI,prctile(PstWLSS(2:end,:),95),'color',col2,'lineWidth',1.1);
plot(cfg.FOI,prctile(PstWLSS(2:end,:),99.9),'color',col1,'lineWidth',1.1);
xlim([cfg.FOI(1) cfg.FOI(end)])
ylabel('Power');
xlabel('Freq. (Hz)');     
set(gca,'TickDir','out');
text(-0.3, 0.15,'Fix. Eff.','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.46,-0.5,'st WLSS','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.44,-1.2,'__________________________','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.2,1.13,'G','units','normalized', 'FontSize', 16);

subplot(6,3,14), hold on; % F.E. max space
plot(cfg.FOI,PstWLSS(1,:),'lineWidth',1.4);
yline(prctile(max(PstWLSS(2:end,:),[],2),95),'color',col1,'lineWidth',1.1,'alpha', 1);
xlim([cfg.FOI(1) cfg.FOI(end)]);   
ylabel('Power');
xlabel('Freq. (Hz)');      
set(gca,'TickDir','out');
text(-0.2,1.13,'H','units','normalized', 'FontSize', 16);

subplot(6,3,15), hold on;%F.E. FDR
pval=sum(PstWLSS(2:end,:)>PstWLSS(1,:),1)./(size(PstWLSS,1)-1);  %p-value  
fdr = calcFDR(pval);
[~,b]=find(fdr<0.1);
plot(cfg.FOI,PstWLSS(1,:),'lineWidth',1.4);
plot(cfg.FOI(b),ones(1,length(b)).*max(PstWLSS(1,:)).*1.1,'*','color',col1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('Power');
xlabel('Freq. (Hz)');  
set(gca,'TickDir','out');
text(-0.2,1.13,'I','units','Normalized', 'FontSize', 16);

subplot(6,3,16), hold on  % R.E. bonferroni
plot(cfg.FOI, tVstWLSS_Ob,'lineWidth',1.4);
plot(cfg.FOI, mean(tVstWLSS_Pr,1),'color',col3,'lineWidth',1.1);
plot(cfg.FOI, prctile(tVstWLSS_Pr,95,1),'color',col2,'lineWidth',1.1);
plot(cfg.FOI, prctile(tVstWLSS_Pr,99.9,1),'color',col1,'lineWidth',1.1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('T-values');
xlabel('Freq. (Hz)');   
set(gca,'TickDir','out');
text(-0.3,0.15,'Rand. Eff.','rotation',90, 'units','normalized', 'FontSize', 11,'FontWeight', 'bold');
text(-0.2,1.13,'J','units','normalized', 'FontSize', 16);

subplot(6,3,17),hold on;  % R.E. max based
plot(cfg.FOI, tVstWLSS_Ob,'lineWidth',1.4);
yline(prctile(max(tVstWLSS_Pr,[],2),95),'color',col1,'lineWidth',1.1,'alpha', 1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('T-values');
xlabel('Freq. (Hz)');     
set(gca,'TickDir','out');
text(-0.2,1.13,'K','units','normalized', 'FontSize', 16);

subplot(6,3,18), hold on,  %F.E. FDR
pval=sum(tVstWLSS_Pr>tVstWLSS_Ob,1)./size(tVstWLSS_Pr,1);  %p-value  
fdr = calcFDR(pval);
[~,b]=find(fdr<0.1);
plot(cfg.FOI,tVstWLSS_Ob,'lineWidth',1.4);
plot(cfg.FOI(b),ones(1,length(b)).*max(tVstWLSS_Ob).*1.1,'*','color',col1);
xlim([cfg.FOI(1) cfg.FOI(end)]);
ylabel('T-values');
xlabel('Freq. (Hz)')  
set(gca,'TickDir','out');
text(-0.2,1.13,'L','units','normalized', 'FontSize', 16);


set(gcf,'Renderer', 'Painters'); 
% orient(gcf,'PORTRAIT')
% print(gcf,'Fig3- Mult Comp Fixed VS Rand Effect','-dpdf')