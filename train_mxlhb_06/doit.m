% This script calls check to check the data and input specifications, transforms the data
% into a more easily useable form, calls hbmcmc to do the iterations, and prints out results.
% Written by Kenneth Train on July 16, 2006, latest edits on Sept 24, 2006.


% Check the input data and specifications

ok=check;
if ok == 1;
    disp('Inputs have been checked and look fine.');
else
    return;
end

% Create Global variables to use in iterations

format short;
SQRTNP=sqrt(NP);

cp=XMAT(:,1); % person

nn=zeros(NCS,1);
for n=1:NCS;
    nn(n,1)=sum(XMAT(:,2) == n,1);
end;
NALTMAX=max(nn);  %Maximum number of alternatives in any choice situation

nn=zeros(NP,1);
for n=1:NP;
   k=(XMAT(:,1)==n);
   k=XMAT(k,2);
   nn(n,1)=1+k(end,1)-k(1,1);
end;
NCSMAX=max(nn);  %Maximum number of choice situations faced by any person

% Data arrays 
% All variables are differenced from the chosen alternative
% Only nonchosen alternatives are included, since V for chosen alt =0
% This reduces number of calculations for logit prob and eliminates need to
% retain dependent variable.

X=zeros(NALTMAX-1,NCSMAX,NV,NP); % Explanatory variables with random coefficients 
%                                  for all choice situations, for each person 
XF=zeros(NALTMAX-1,NCSMAX,NF,NP); % Explanatory variables with fixed coefficients 
%                                  for all choice situations, for each person 
S=zeros(NALTMAX-1,NCSMAX,NP); % Identification of the alternatives in each choice situation, for each person

for n=1:NP;  %loop over people
 cs=XMAT(cp == n,2);
 yy=XMAT(cp == n,3);
 if NV > 0
    xx=XMAT(cp == n, IDV(:,1));
 end
 if NF > 0
    xxf=XMAT(cp == n, IDF(:,1));
 end
 t1=cs(1,1);
 t2=cs(end,1);
 for t=t1:t2; %loop over choice situations
     k=sum(cs==t)-1; %One less than number of alts = number of nonchosen alts
     S(1:k,1+t-t1,n)=ones(k,1);
     if NV>0
        X(1:k,1+t-t1,:,n)=xx(cs==t & yy == 0,:)-repmat(xx(cs==t & yy == 1,:),k,1);
     end
     if NF>0
        XF(1:k,1+t-t1,:,n)=xxf(cs==t & yy == 0,:)-repmat(xxf(cs==t & yy == 1,:),k,1);
     end;
 end
end

%S=sparse(S);
clear XMAT cp cs yy xx t1 t2 k
    
randn('state',SEED1)  %Though nothing in the code is random yet
rand('state',SEED1)   %For draws from uniform


[adraws,ddraws,cmeans,fdraws]=hbmcmc;

%Save if user specifies
if KEEPMN == 1
    if NV > 0
      save(PUTA,'adraws');
      save(PUTD,'ddraws');
    end
    if NF > 0
      save(PUTF,'fdraws');
    end
end
if WANTINDC == 1 & NV>0
    save(PUTC,'cmeans');
end
disp('RESULTS');
disp(' ');
disp(' ')
if NF>0
disp('FIXED COEFFICIENTS');
disp('F-hat, se(F-hat)')
disp('From a classical perspective these are estimates and standard errors. ');
disp('From a Bayesian perspective these are mean and std dev of posterior. ');
disp(' ');
disp(' ');
fhat=mean(fdraws,2);
fsd=std(fdraws,0,2);
disp('                      F      ');
disp('              ------------------ ');
disp('              Est/Mean   SE/StDv ');
for r=1:length(NAMESF);
    fprintf('%-10s %10.4f %10.4f\n', NAMESF{r,1}, [fhat(r,1) fsd(r,1)]);
end
disp(' ');
end;

if NV>0;
disp('RANDOM COEFFICIENTS');
disp('A-hat, se(A-hat), diag(D-hat), se(diag(D-hat).');
disp('From a classical perspective these are estimates and standard errors. ');
disp('From a Bayesian perspective these are mean and std dev of posterior. ');
disp(' ');
disp(' ');

ahat=mean(adraws,2);
dhat=mean(ddraws,3);
asd=std(adraws,0,2);
dsd=std(ddraws,0,3);
disp('                      A            Diagonal elements of D');
disp('              ------------------   -----------------------');
disp('              Est/Mean   SE/StDv    Est/Mean   SE/StDv');
for r=1:length(NAMES);
    fprintf('%-10s %10.4f %10.4f %10.4f %10.4f\n', NAMES{r,1}, [ahat(r,1) asd(r,1) dhat(r,r) dsd(r,r)]);
end
disp(' ');
if FULLCV == 1
    disp('D-hat: Estimate of D or mean of posterior of D:');
    disp(dhat);
    disp('Correlations in D-hat');
    disp((dhat./repmat(sqrt(diag(dhat)),1,NV))./repmat(sqrt(diag(dhat))',NV,1));
end

%Create draws of coefficients from A-hat and D-hat
randn('state',SEED2);
cdhat=chol(dhat);
B=repmat(ahat,1,NRFOR)+cdhat'*randn(NV,NRFOR);
C=trans(B);
disp('Distribution of coefficients in population implied by A-hat and D-hat.');
disp(['Uses ' num2str(NRFOR)  ' draws from N(A-hat,D-hat) and transforms each of them as']);
disp('specified in IDV to create draws of coefficients.');
disp(' ');
jj={'normal';'lognormal';'truncnormal';'S_B';'normal0mn'};
disp('                            Mean      StdDev     Share<0    Share=0');
kk=[mean(C,2) std(C,0,2) mean((C < 0),2) mean((C == 0),2)];
for r=1:length(NAMES);
    mm=IDV(r,2);
    fprintf('%-10s %-11s %10.4f %10.4f %10.4f %10.4f\n', NAMES{r,1}, jj{mm,1},kk(r,:));
end
disp(' ');
if FULLCV == 1;
    disp('Correlation among coefficients');
    disp(corrcoef(C'));
    disp(' ');
end
end
%Calculate simulated log-likelihood at A-hat, D-hat, F-hat

p=zeros(1,NP);
for r=1:NRLL
    if NV>0
        B=repmat(ahat,1,NP)+cdhat'*randn(NV,NP);
    else
        B=zeros(NV,NP);
    end
    if NF == 0
        fhat=zeros(NF,1);
    end
    p=p+logit(fhat,trans(B)); 
end
p=p./NRLL;
disp('Simulated log-likelihood value at A-hat, D-hat, and F-hat');
disp(['based on ' num2str(NRLL) ' draws per person for the random coefficients:']);
disp(sum(log(p),2));



    

