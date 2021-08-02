%% Simple code to fit QR, Pss and GP models, all the parameters assumed to vary linearly in time
% Mdl(xi, sgm, rho, psi) with xi=xi_0 + Tim * xi_1 etc
% Will run as is for toy data to check 
% Need to input data as structure (see occurrences of USER INPUT) below
%
% QR = non-stationary Quantile Regression threshold
% Pss = non-stationary PoiSSon could model for threshold exceedances
% GP = Generalised Pareto model for size of threshold exceedances

%% Output for further plotting / investigation etc
% Output from the analysis is saved in a structure C within file MCMC.mat
% A typical structure C (when 8 different threshold non-exceedance probabilities are used) is
%
%% C
% 
%       Nep: [8×1 double] Non-exceedance probabilities for thresholds 
%      nNep: 8            Number of NEPs
%        QR: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}  QR model details (inc posterior sample)
%       Pss: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}  Pss model details (inc posterior sample)
%        GP: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}  GP model details (inc posterior sample)
%        RV: {[1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]  [1×1 struct]}  RV details (inc posterior sample)
%     RVCmp: [1×1 struct]                                                                                                      RV comparison summary
%    PrmSmm: [1x1 struct]                                                                                                      Assessment of slope parameter changes
%
%% C.QR{q}, C.Pss{q}, C.GP{q} for q=1,2,..., nNep are structures like
% 
%        Lkl: 'QR'
%       nPrm: 2
%     PrmNms: {2×1 cell}
%        Nep: 0.6000
%       nItr: 10000
%      n2Plt: 5000
%     NgtStr: 0.1000
%     AdpItr: 1000
%     AdpBet: 0.0500
%     PrmStr: [2×1 double]
%     AccRat: [10000×2 double]
%        Prm: [10000×2 double]
%        Nll: [10000×1 double]
%     PrmUse: [2×1 double]
%
%% Key output for further plotting etc are
%
% C.QR{q}.Prm   nItr x 2 values of psi0 and psi1 from MCMC (in general it is safe to use the last 9000; first 1000 might involve "burn-in")
% C.Pss{q}.Prm  nItr x 2 values of rho0 and rho1 from MCMC (in general it is safe to use the last 9000; first 1000 might involve "burn-in")
% C.GP{q}.Prm   nItr x 4 values of xi0, xi1, sigma0 and sigma1 from MCMC (in general it is safe to use the last 9000; first 1000 might involve "burn-in")
% C.RV{q}.Prm   nRls x 2 values of return values at "time zero" and "time end" generated using C.QR.Prm, C.Pss.Prm and C.GP.Prm
% C.RVCmp.Prb   nNep x 1 probabilities that the 100-year maximum at "time end" is larger than at "time zero" for each NEP
% C.RVCmp.Qnt   nNep x 6 Quantlies (2.5%, 50% and 97.5%) at "time end" and quantlies (2.5%, 50% and 97.5%) at "time end" for each NEP
% C.PrmCmp.Prb  nNep x 4 probabilities that slope terms psi1, rho1, xi1 and sigma1 are "significant" for each NEP
% C.PrmCmp.Qnt  nNep x 3 x 4 Quantlies (2.5%, 50% and 97.5%) for psi1, rho1, xi1 and sigma1 at each NEP
% C.txt         A text file summarising information in C.RVCmp.Prb and C.PrmSmm.Prb

%% Set up
% clc; clear; clf;
VrbNms={'$\xi$';'$\sigma$';'$\psi$'};

%% Simulate a sample of data
if 1; %for testing

    % True parameters P0=[xi0;xi1;sgm0;sgm1;rho0;rho1;psi0;psi1] of linear regression 
    
    %% *** USER INPUT *** Pick the type of simulated data
    %Tst='Up';   %Upward trend
    Tst='Same';  %No trend
    %Tst='Down'; %Downward trend
    
    switch Tst
        case 'Up';
            X.Prm0=[...
                -0.3;0.05;
                2;0.1;...
                20;1;...
                2;0;...
                ];
        case 'Same';
            X.Prm0=[...
                -0.3;0;
                1;0;...
                20;0;...
                2;0;...
                ];
        case 'Down';
            X.Prm0=[...
                -0.3;0;
                3;-0.5;...
                20;0;...
                2;0;...
                ];
    end;
    
    % Time variable
    % NB A COMMON time value used for observations in the same year
    X.nYr=85;
    tYr=(1:X.nYr)';  % Time in years

    % True parameter estimates per year
    X.XSM0=[ones(X.nYr,1)*X.Prm0(1)+(tYr/X.nYr)*X.Prm0(2) ones(X.nYr,1)*X.Prm0(3)+(tYr/X.nYr)*X.Prm0(4) ...
        ones(X.nYr,1)*X.Prm0(5)+(tYr/X.nYr)*X.Prm0(6) ones(X.nYr,1)*X.Prm0(7)+(tYr/X.nYr)*X.Prm0(8)];

    % Number of occurrences per annum
    tOcc=poissrnd(X.XSM0(:,3));
    
    % Generate data from GP
    k=0;
    X.nT=sum(tOcc);
    X.Tim=nan(X.nT,1);
    X.Dat=nan(X.nT,1);
    for iY=1:X.nYr;
        for iO=1:tOcc(iY);
            k=k+1;
            X.Tim(k)=tYr(iY)/X.nYr;
            X.Dat(k)=gprnd(X.XSM0(iY,1),X.XSM0(iY,2),X.XSM0(iY,4));
        end;
    end;
    
    X, % See the structure
    
    subplot(2,2,1); plot(tYr,tOcc,'ko');
    subplot(2,2,2); plot(X.Tim*X.nYr,X.Dat,'ko');
    
end;

%% ***USER INPUT*** Read in your data here
if 0;
    %X.nYr ; %  1   x 1 number of years
    %X.nT  ; %  1   x 1 number of occurrences
    %X.Tim ; % nT   x 1 years on [0,1] (so that floor((X.Tim*X.nYr)+1) gives the year number
    %X.Dat ; % nT   x 1 data
    load userInput.mat;
end;

%% ***USER INPUT*** Specify NEPs to consider
if 1;
    %C.Nep=(0.6:0.05:0.95)'; % (0.6:0.05:0.9)' is a good range; but maybe you want to use (0.7:0.1:0.9)' to get going
    %C.Nep=(0.7:0.1:0.9)';
    C.Nep=[0.7;0.8];
    C.nNep=size(C.Nep,1);
end;

%% Estimate extreme value threshold (linear Quantile Regression)
if 1;
    
    for iN=1:C.nNep
        
        C.QR{iN}.Lkl='QR';       % Likelihood
        C.QR{iN}.nPrm=2;         % Number of parameters
        C.QR{iN}.PrmNms={'$\psi_0$';'$\psi_1$';}; % Names for parameters
        C.QR{iN}.Nep=C.Nep(iN);  % NEP

        C.QR{iN}.nItr=10000;     % Number of MCMC iterations - 1e4 minipsim when used in anger
        C.QR{iN}.n2Plt=5000;     % Number of iterations from end of chain to "beleive"
        
        C.QR{iN}.NgtStr=0.1;     % Candidate random walk standard deviation - don't change
        C.QR{iN}.AdpItr=1000;    % Number of warm up iterations - don't change
        C.QR{iN}.AdpBet=0.05;    % Adaptive MC - don't change        C.Nep=X.Nep(iN);
        
        C.QR{iN}.PrmStr=[quantile(X.Dat,C.Nep(iN));0]; % Constant starting solution for quantile regression
        
        C.QR{iN}=Mcmc(X,C.QR{iN});   % Run MCMC algorithm
        
        C.QR{iN}.PrmUse=mean(C.QR{iN}.Prm(C.QR{iN}.nItr-C.QR{iN}.n2Plt+1:C.QR{iN}.nItr,:))'; % Use posterior mean for subsequent inference
        
        tStr=sprintf('GpNonStt-Mdl%s-Nep%g',C.QR{iN}.Lkl,C.QR{iN}.Nep); pDatStm(tStr); pGI(tStr,2); % Save plot  
        tFil=sprintf('MCMC'); save(tFil,'C'); % Save whole chain

    end;
    
end;

%% Estimate rate of threshold exceedance per annum (linear Poisson Process)
if 1;
        
    for iN=1:C.nNep
        
        C.Pss{iN}.Lkl='Pss';      % Likelihood
        C.Pss{iN}.Nep=C.Nep(iN);  % NEP
        C.Pss{iN}.PrmNms={'$\rho_0$';'$\rho_1$';}; % Names for parameters
        C.Pss{iN}.nPrm=2;         % Number of parameters
        
        C.Pss{iN}.nItr=10000;     % Number of MCMC iterations - 1e4 minipsim when used in anger
        C.Pss{iN}.n2Plt=5000;     % Number of iterations from end of chain to "beleive"
        
        C.Pss{iN}.NgtStr=0.1;     % Candidate random walk standard deviation - don't change
        C.Pss{iN}.AdpItr=1000;    % Number of warm up iterations - don't change
        C.Pss{iN}.AdpBet=0.05;    % Adaptive MC - don't change        C.Nep=X.Nep(iN);
        
        % Estimate Poisson count for threshold exceedances
        t1=(X.Dat-(ones(X.nT,1)*C.QR{iN}.PrmUse(1)+X.Tim*C.QR{iN}.PrmUse(2)))>0; % Threshold exceedances
        for iY=1:X.nYr;
            t2=floor(X.Tim*X.nYr)>=(iY-1) & floor(X.Tim*X.nYr)<iY; % Particular year
            C.Pss{iN}.Cnt(iY,:)=sum(t1(t2==1));
            C.Pss{iN}.CntTim(iY,:)=(iY-0.5)/X.nYr; % Take middle of year
        end;
        
        % Constant starting solution from Poisson fit
        C.Pss{iN}.PrmStr=[poissfit(C.Pss{iN}.Cnt);0]; 
        
        C.Pss{iN}=Mcmc(X,C.Pss{iN});   % Run MCMC algorithm
        
        tStr=sprintf('GpNonStt-Mdl%s-Nep%g',C.Pss{iN}.Lkl,C.Pss{iN}.Nep); pDatStm(tStr); pGI(tStr,2); % Save plot
        tFil=sprintf('MCMC'); save(tFil,'C'); % Save whole chain

    end;
    
end;

%% Estimate size of threshold exceedance per annum (linear Generalised Pareto)
if 1;
            
    for iN=1:C.nNep
        
        C.GP{iN}.Lkl='GP'; 
        C.GP{iN}.Nep=C.Nep(iN);  % NEP
        C.GP{iN}.PrmNms={'$\xi_0$';'$\xi_1$';'$\sigma_0$';'$\sigma_1$';}; % Names for parameters
        C.GP{iN}.nPrm=4;     % Number of parameters
        
        C.GP{iN}.nItr=10000;     % Number of MCMC iterations - 1e4 minipsim when used in anger
        C.GP{iN}.n2Plt=5000;     % Number of iterations from end of chain to "beleive"
        
        C.GP{iN}.NgtStr=0.1;     % Candidate random walk standard deviation - don't change
        C.GP{iN}.AdpItr=1000;    % Number of warm up iterations - don't change
        C.GP{iN}.AdpBet=0.05;    % Adaptive MC - don't change        C.Nep=X.Nep(iN);
        
        % Isolate threshold exceedances and times of occurrence
        t1=X.Dat-(ones(X.nT,1)*C.QR{iN}.PrmUse(1)+X.Tim*C.QR{iN}.PrmUse(2)); % Threshold exceedances
        C.GP{iN}.Exc=t1(t1>0);
        C.GP{iN}.ExcTim=X.Tim(t1>0);
        
        % Constant starting solution from GP fit
        t=gpfit(C.GP{iN}.Exc); 
        if t(1)<-0.5; t(1)=-0.4; end;
        if t(1)>0.5; t(1)=0.4; end;
        C.GP{iN}.PrmStr=[t(1);0;t(2);0]; 
        
        C.GP{iN}=Mcmc(X,C.GP{iN});   % Run MCMC algorithm
        
        tStr=sprintf('GpNonStt-Mdl%s-Nep%g',C.GP{iN}.Lkl,C.GP{iN}.Nep); pDatStm(tStr); pGI(tStr,2); % Save plot  
        tFil=sprintf('MCMC'); save(tFil,'C'); % Save whole chain

    end;
    
end;

%% Compare 100-year values at start and end of study as a function of threshold
if 1;

    %% Calculate return values
    for iN=1:C.nNep
     
        C.RV{iN}.RtrPrd=100; % Return period of interest
        C.RV{iN}.nRls=1000;  % Number of realisations to use (1000 is good)
        
        % Parameter estimates at start
        t1=randi(C.QR{iN}.nItr-C.QR{iN}.n2Plt,C.RV{iN}.nRls,1)+C.QR{iN}.n2Plt;
        t2=randi(C.Pss{iN}.nItr-C.Pss{iN}.n2Plt,C.RV{iN}.nRls,1)+C.Pss{iN}.n2Plt;
        t3=randi(C.GP{iN}.nItr-C.GP{iN}.n2Plt,C.RV{iN}.nRls,1)+C.GP{iN}.n2Plt;
        if 0; %This uses the MEAN threshold (and so does not propagate threshold uncertainty) ***USER INPUT***
            tPsi=C.QR{iN}.PrmUse(1);
        else; %This uses a random threshold (and so "slightly" overestimates the effect of threshold uncertainty) ***USER INPUT***
            tPsi=C.QR{iN}.Prm(t1,1);
        end;
        tRat=C.Pss{iN}.Prm(t2,1)+C.Pss{iN}.CntTim(1)*C.Pss{iN}.Prm(t2,2); % Rate is per annum, modelled at midpoint of year
        tXi=C.GP{iN}.Prm(t3,1);
        tSgm=C.GP{iN}.Prm(t3,3);      
        C.RV{iN}.RV(:,1)=(tSgm./tXi).*( (C.RV{iN}.RtrPrd.*tRat).^tXi - 1) + tPsi; % Return value for year zero
        
        % Parameter estimates at end
        if 0; %0 hear means that the SAME random shots are used to compare end with start; this may be advantageous ***USER INPUT***
              %In pGevNonStt analysis, we use the SAME random shot to compare end with start => we should have 0 here for fair comparison
            t1=randi(C.QR{iN}.nItr-C.QR{iN}.n2Plt,C.RV{iN}.nRls,1)+C.QR{iN}.n2Plt;
            t2=randi(C.Pss{iN}.nItr-C.Pss{iN}.n2Plt,C.RV{iN}.nRls,1)+C.Pss{iN}.n2Plt;
            t3=randi(C.GP{iN}.nItr-C.GP{iN}.n2Plt,C.RV{iN}.nRls,1)+C.GP{iN}.n2Plt;
        end;
        if 0; %This uses the MEAN threshold (and so does not propagate threshold uncertainty) ***USER INPUT***
            tPsi=sum(C.QR{iN}.PrmUse,2);
        else; %This uses a random threshold (and so "slightly" overestimates the effect of threshold uncertainty) ***USER INPUT***
            tPsi=sum(C.QR{iN}.Prm(t1,:),2);
        end;
        tRat=C.Pss{iN}.Prm(t2,1)+C.Pss{iN}.CntTim(end)*C.Pss{iN}.Prm(t2,2); % Rate is per annum, modelled at midpoint of year
        tXi=sum(C.GP{iN}.Prm(t3,1:2),2);
        tSgm=sum(C.GP{iN}.Prm(t3,3:4),2);      
        C.RV{iN}.RV(:,2)=(tSgm./tXi).*( (C.RV{iN}.RtrPrd.*tRat).^tXi - 1) + tPsi; % Return value for last year
        
        % Summary statistics of differences
        C.RVCmp.Prb(iN,:)=mean(C.RV{iN}.RV(:,2)>C.RV{iN}.RV(:,1)); % Prob. that RVEnd > RVStart
        C.RVCmp.Qnt(iN,:)=[quantile(C.RV{iN}.RV(:,1),[0.025 0.5 0.975]) quantile(C.RV{iN}.RV(:,2),[0.025 0.5 0.975])]; % Quantiles
        
    end;
    
    %% Figure
    clf;
    subplot(1,2,1); 
    plot(C.Nep,C.RVCmp.Prb,'ko-'); 
    pAxsLmt; pDflHug;
    title 'Prb(RVEnd$>$RVStart)'; xlabel 'Threshold NEP';
    subplot(1,2,2); hold on;
    plot(C.Nep,C.RVCmp.Qnt(:,1),'ko-'); 
    plot(C.Nep,C.RVCmp.Qnt(:,2),'ko-','linewidth',3); 
    plot(C.Nep,C.RVCmp.Qnt(:,3),'ko-'); 
    plot(C.Nep,C.RVCmp.Qnt(:,4),'ro-'); 
    plot(C.Nep,C.RVCmp.Qnt(:,5),'ro-','linewidth',3); 
    plot(C.Nep,C.RVCmp.Qnt(:,6),'ro-'); 
    if isfield(X,'Prm0')==1; % True values are known
        tPsi=X.Prm0(7);
        tRat=X.Prm0(5)+C.Pss{1}.CntTim(1)*X.Prm0(6);
        tSgm=X.Prm0(3);
        tXi=X.Prm0(1);
        tRVTru=(tSgm./tXi).*( (C.RV{1}.RtrPrd.*tRat).^tXi - 1) + tPsi;
        plot(C.Nep,ones(C.nNep,1)*tRVTru,'k--');
        tPsi=sum(X.Prm0(7:8));
        tRat=X.Prm0(5)+C.Pss{1}.CntTim(end)*X.Prm0(6);
        tSgm=sum(X.Prm0(3:4));
        tXi=sum(X.Prm0(1:2));
        tRVTru=(tSgm./tXi).*( (C.RV{1}.RtrPrd.*tRat).^tXi - 1) + tPsi;
        plot(C.Nep,ones(C.nNep,1)*tRVTru,'r--');
    end;
    pAxsLmt; pDflHug;
    title 'Quantiles of RVStart (k), RVEnd (r) [True ---]'; xlabel 'Threshold NEP';
    pAxsLmt; pDflHug;
    pDatStm; pGI('GpNonStt-100YearReturnValue',2);
    
    %% Summary statistics to screen
    clc;
    fid=fopen('C.txt','w+');
    fprintf(1,'SUMMARY FOR CHANGE IN RETURN VALUES\n');
    for iN=1:C.nNep;
        if C.RVCmp.Prb(iN)>0.975 || C.RVCmp.Prb(iN)<=0.025;
            fprintf(1,'NEP=%4.2f Prb(RVEnd>RVStart)=%4.2f SIGNIFICANT\n',C.Nep(iN),C.RVCmp.Prb(iN));
            fprintf(fid,'NEP=%4.2f Prb(RVEnd>RVStart)=%4.2f SIGNIFICANT\n',C.Nep(iN),C.RVCmp.Prb(iN));
        else;
            fprintf(1,'NEP=%4.2f Prb(RVEnd>RVStart)=%4.2f (not significant)\n',C.Nep(iN),C.RVCmp.Prb(iN));
            fprintf(fid,'NEP=%4.2f Prb(RVEnd>RVStart)=%4.2f (not significant)\n',C.Nep(iN),C.RVCmp.Prb(iN));
        end;
    end;    
    fclose(fid);
    
end;

%% Plot parameter estimates nicely
if 1;
    
    clf;
    tTtl={'QR';'Pss';'GP shape';'GP scale'};
    for iN=1:C.nNep;
        for j=1:8;
            switch j
                case 1; tDat=C.QR{iN}.Prm(C.QR{iN}.n2Plt+1:end,1); tTxt=C.QR{iN}.PrmNms{1};
                case 2; tDat=C.Pss{iN}.Prm(C.Pss{iN}.n2Plt+1:end,1); tTxt=C.Pss{iN}.PrmNms{1};
                case 3; tDat=C.GP{iN}.Prm(C.GP{iN}.n2Plt+1:end,1); tTxt=C.GP{iN}.PrmNms{1};
                case 4; tDat=C.GP{iN}.Prm(C.GP{iN}.n2Plt+1:end,3); tTxt=C.GP{iN}.PrmNms{3};
                case 5; tDat=C.QR{iN}.Prm(C.QR{iN}.n2Plt+1:end,2); tTxt=C.QR{iN}.PrmNms{2};
                case 6; tDat=C.Pss{iN}.Prm(C.Pss{iN}.n2Plt+1:end,2); tTxt=C.Pss{iN}.PrmNms{2};
                case 7; tDat=C.GP{iN}.Prm(C.GP{iN}.n2Plt+1:end,2); tTxt=C.GP{iN}.PrmNms{2};
                case 8; tDat=C.GP{iN}.Prm(C.GP{iN}.n2Plt+1:end,4); tTxt=C.GP{iN}.PrmNms{4};
            end;
            subplot(2,4,j); hold on;
            pDns(tDat,10,pClr(iN));
            if iN==C.nNep
                if j>4;
                    t=get(gca,'ylim');
                    plot(ones(2,1)*0, [0;t(2)*0.9],'color','k','linewidth',3);
                end;
                if j==1;
                    ylabel('Intercept terms');
                elseif j==5;
                    ylabel('Slope terms');
                end;
                if j<=4;
                    title(tTtl{j});
                end;
                xlabel(tTxt,'interpreter','LaTeX');
                pAxsLmt; pDflHug;
            end;
        end;    
    end;
    tClr={'r';'m';'g';'c';'b';'k';'gr';'br'};
    t='NEPs: ';
    for iN=1:C.nNep;
        t=sprintf('%s %g (%s)',t,C.Nep(iN),tClr{iN});
    end;
    pDatStm(t);
    pGI('GpNonStt-Parameters',2);
       
    %% Summary statistics to screen
    fid=fopen('C.txt','w+');
    fprintf(1,'\nSUMMARY FOR SIGNIFICANT SLOPE PARAMETERS\n');
    for j=1:4;
        for iN=1:C.nNep;
            switch j
                case 1; tDat=C.QR{iN}.Prm(C.QR{iN}.n2Plt+1:end,2); tTxt=C.QR{iN}.PrmNms{2};
                case 2; tDat=C.Pss{iN}.Prm(C.Pss{iN}.n2Plt+1:end,2); tTxt=C.Pss{iN}.PrmNms{2};
                case 3; tDat=C.GP{iN}.Prm(C.GP{iN}.n2Plt+1:end,2); tTxt=C.GP{iN}.PrmNms{2};
                case 4; tDat=C.GP{iN}.Prm(C.GP{iN}.n2Plt+1:end,4); tTxt=C.GP{iN}.PrmNms{4};
            end;
            tPrb=sum(tDat>0)/size(tDat,1);
            tSum=min(sum(tDat>0),sum(tDat<0))/size(tDat,1);        
            tQnt=quantile(tDat,[0.025 0.5 0.975]);
            if tSum <0.025;
                fprintf(1,'%s: NEP=%4.2f Prb(Prm>0)=%4.2f SIGNIFICANT\n',tTxt, C.Nep(iN),tPrb);
                fprintf(fid,'%s: NEP=%4.2f Prb(Prm>0)=%4.2f SIGNIFICANT\n',tTxt, C.Nep(iN),tPrb);
            else;
                fprintf(1,'%s: NEP=%4.2f Prb(Prm>0)=%4.2f (not significant)\n',tTxt, C.Nep(iN),tPrb);
                fprintf(fid,'%s: NEP=%4.2f Prb(Prm>0)=%4.2f (not significant)\n',tTxt, C.Nep(iN),tPrb);
            end;
            C.PrmSmm.Prb(iN,j)=tPrb;
            C.PrmSmm.Qnt(iN,:,j)=tQnt;
        end;
        fprintf(1,'\n');
        fprintf(fid,'\n');
    end;
    fclose(fid);
    
end;

%% Update output file
if 1;
    tFil=sprintf('MCMC'); save(tFil,'C'); % Save whole chain
end;