function [ll,g]=flexll(alpha)

global PROBS Z NZ NDRAWS IDV IDNR XMAT BETAS NNR NP NCS NROWS

% PROBS: NPxNDRAWS matrix of probabilities for each person at each sampled coefficients
% Z:     NPxNDRAWSxNZ array of variables that explain distribution of coefficients
% NZ: scalar number of variables in Z
% NDRAWS: scalar number of draws in PROBS

% Input alpha is NZx1 vector of coefficients

v_ran=sum(bsxfun(@times,XMAT(:,IDV),BETAS(XMAT(:,1),:,:)),2); 
CoefNR = reshape(alpha(end-NNR+1:end), 1, NNR);
BETAS_FIX=repmat(CoefNR,[NP,1,NDRAWS]); %NPxNNRxNDRAWS
v_fix=sum(bsxfun(@times,XMAT(:,IDNR),BETAS_FIX(XMAT(:,1),:,:)),2); 
v = v_ran + v_fix;
v=squeeze(v);  %NROWSxNDRAWS
v=exp(v);
sparsematrix=bsxfun(@eq,sparse(1:NCS)',XMAT(:,2)'); %NCSxNROWS
denom=double(sparsematrix)*v;                           %NCSxNDRAWS
p=v(XMAT(:,3)==1,:)./denom;                     %NCSxNDRAWS
personid=XMAT(XMAT(:,3)==1,1);                  %NCSx1
sparsematrix=bsxfun(@eq,sparse(1:NP)',personid');    %NPxNCS
PROBS=double(sparsematrix)*log(p);                      %NPxNDRAWS
PROBS=exp(PROBS);                               %NPxNDRAWS

w=sum(bsxfun(@times,Z,reshape(alpha(1:NZ),[1,1,NZ])),3);  %NPxNDRAWS
w(w<-500)=-500;  %As precaution against extreme parameters
w(w>500)=500;
w=squeeze(exp(w));
w=w./repmat(sum(w,2),1,NDRAWS); %NPxNDRAWS
logit_w=w.*PROBS; %NPxNDRAWS
mix_probs=sum(logit_w,2); %NPx1
ll=sum(log(mix_probs),1); %1x1
ll=-ll; %To minimize


if nargout>1;
  condw=bsxfun(@times,logit_w,1./mix_probs); %NPxNDRAWS
  g1=bsxfun(@times,(condw-w),Z); %NPxNDRAWSxNZ
  g1=sum(sum(g1,1),2); %1x1xNZ
  g1=reshape(g1,NZ,1); %NZx1
  
  NRX =XMAT(XMAT(:,3)==1,IDNR);
  sparsematrix=bsxfun(@eq,sparse(1:NP)',personid');    %NPxNCS
  T1 = double(sparsematrix)*NRX;
  T1 =repmat(reshape(T1,[NP,1,NNR]),[1,NDRAWS,1]);
  
  v_new =repmat(reshape(v,[NROWS,1,NDRAWS]),[1,NNR,1]);
  X_V = bsxfun(@times,v_new,XMAT(:,IDNR)); %NROWSxNNRxNDRAWS
  sparsematrix=bsxfun(@eq,sparse(1:NCS)',XMAT(:,2)'); %NCSxNROWS
  
  numerator = zeros(NCS,NDRAWS,NNR);
  for i=1:NNR
      numerator(:,:,i) = double(sparsematrix)*squeeze(X_V(:,i,:));  
  end        
  ratio = bsxfun(@times,numerator,1./denom); % NCSxNDRAWSxNNR
  sparsematrix=bsxfun(@eq,sparse(1:NP)',personid');    %NPxNCS
  
  T2 = zeros(NP,NDRAWS,NNR);
  for i=1:NNR
      T2(:,:,i) = double(sparsematrix)*ratio(:,:,i);  
  end   
  ZNR = T1-T2;
  g2=bsxfun(@times,condw,ZNR); %NPxNDRAWSxNNR
  g2=sum(sum(g2,1),2); %1x1xNNR
  g2=reshape(g2,NNR,1); %NNRx1
   
  g=-[g1;g2]; %To minimize
end;
  
  
