function C=Mcmc(Y,C);
% function C=Mcmc(Y,C);
%
% MCMC for GEV fit, all parameters take form p_0 + time * p_1
% Objective is to estimate xi0, xi1, sigma0, sigma1, mu0 and mu1

%% Parameter names


%% Loop over MCMC iterations
C.AccRat=nan(C.nItr,C.nPrm);
for iI=1:C.nItr;
    
    %% Find valid starting solution if iI=1
    if iI==1;
        
        NllStr=EstNll(C.PrmStr,Y,C);
        
        if isnan(NllStr)==1; %no valid starting solution found. Terminate.
            fprintf(1,'Warning: invalid starting solution. Terminating.\n');
            return;
        else %make the current state the starting state for MCMC
            Prm=C.PrmStr;
            Nll=NllStr;
            fprintf(1,'Starting solution found. Starting MCMC for case %s with NEP %g\n',C.Lkl,C.Nep);
        end
        
    end;
    
    %% Iteration counter on screen
    if rem(iI,100)==0;
        fprintf(1,'+');
    elseif rem(iI,10)==0;
        fprintf(1,'.');
    end;
    if rem(iI,1000)==0;
        fprintf(1,'\n');
    end;
    
    %% Loop over parametric forms
    for iP=1:C.nPrm; %Metropolis Hastings in Gibbs, one parameter at a time
        
        %% Define candidate in terms of current state
        PrmC=Prm;
        
        if iI<=C.AdpItr; %fixed nugget
            PrmC(iP)=PrmC(iP)+randn*C.NgtStr;
        elseif iP==1; %adaptive Metropolis for A B M S
            jP=[1:C.nPrm]'; nJ=size(jP,1);
            t1=real((1-C.AdpBet)*2.38*sqrtm(cov(C.Prm(max(1,iI-999):iI-1,jP))/nJ)*randn(nJ,1));
            t2=C.AdpBet*0.1*(randn(nJ,1)/nJ);
            PrmC=PrmC+t1+t2;
            C.AccRat(iI-1,jP(2:end))=NaN;
        end;
        
        %% Evaluate likelihood at current state and candidate
        NllC=EstNll(PrmC,Y,C);
        
        if isreal(NllC)==0;
            fprintf(1,'Imaginary NllC\n');
        end
        
        %% MH acceptance step
        if (exp(-NllC+Nll) > rand) && isinf(NllC)==0 && isnan(NllC)==0;
            Prm=PrmC;
            Nll=NllC;
            if iI==1;
                C.AccRat(iI,iP)=1;
            else;
                if iI>100; %Only use last 100 iterations to adjust acceptance rate
                    jI=100;
                else;
                    jI=iI;
                end;
                if iI<=C.AdpItr;
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                elseif iP==1;
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1)+1)/jI;
                else;
                    C.AccRat(iI,iP)=NaN;
                end;
            end;
        else;
            if iI==1; %Only use last 100 iterations to adjust acceptance rate
                C.AccRat(iI,iP)=0;
            else;
                if iI>100; %Only use last 100 iterations to adjust acceptance rate
                    jI=100;
                else;
                    jI=iI;
                end;
                if iI<=C.AdpItr;
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                elseif iP==1;
                    C.AccRat(iI,iP)=(C.AccRat(iI-1,iP)*(jI-1))/jI;
                else;
                    C.AccRat(iI,iP)=NaN;
                end;
            end;
        end;
        
    end;
    
    % Update after complete iteration over variables
    C.Prm(iI,:)=Prm';
    C.Nll(iI,:)=Nll;
    
    
    if rem(iI,1000)==0;
        
        %% Trace plots
        if C.nPrm==2;
            tc=1;
        elseif C.nPrm==4;
            tc=2;
        else;
            tc=ceil(C.nPrm/2);
        end;
        for j=1:C.nPrm;
            subplot(2,tc+1,j);
            plot(C.Prm(:,j),'k-');
            title(C.PrmNms{j},'interpreter','latex');
            pAxsLmt; pDflHug;
        end;
        subplot(2,tc+1,2*tc+1); plot(C.Nll,'k-'); title 'NLL'; pAxsLmt; pDflBig;
        subplot(2,tc+1,2*tc+2); plot(C.AccRat); title 'Acceptance rates'; pAxsLmt; pDflBig;
        drawnow;
                
    end;
        
end;

return;