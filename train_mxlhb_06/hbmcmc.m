% This function performs the iteration of draws of A,D,and B.
% It prints iteration results, keeps up with when burn-in is over,
% and retains draws as specified after burn-in.
% Written by Kenneth Train on July 17, 2006.

%It returns three outputs:
%   adraws contains the draws of A. It is a matrix of dimension NVxNEREP where
%      each column is a draw of vector A.
%   ddraws contains the draws of D. It is an array of dimension NVxNVxNEREP where
%      each ddraws(:,:,r) is a draw of matrix D.
%   cmeans contains the means of the draws of C (which are transformations of the draws of B). 
%      It is a matrix of dimension NVxMP.


function [adraws,ddraws,cmeans,fdraws]=hbmcmc

global A B C D P NV NP
global RHO FULLCV NCREP INFOSKIP NEREP NSKIP
global F NF RHOF

adraws=zeros(NV,NEREP);
ddraws=zeros(NV,NV,NEREP);
fdraws=zeros(NF,NEREP);
cmeans=zeros(NV,NP);
P=logit(F,trans(B));

% Iterate 

kept=0;
for r=1:(NCREP+(NEREP.*NSKIP))
   if NV>0   %For random coefficients
     [B,P,RHO,accept]=nextB;
     A=nextA;
     if FULLCV == 1
        D=nextD;
     else
        D=nextDdiag;
     end
   end
   if NF>0 %For fixed coefficients
       [F,P,keep]=nextF;
       kept=kept+keep;
   end
   
   if ceil(r./INFOSKIP) == (r./INFOSKIP)
      disp('Iteration' );
      disp(r);
      if NV>0
          disp('Random coefficients');
          disp('RHO       Accept rate');
          disp([RHO accept]);
          disp('A and Means of B and C:');
          disp([A mean(B,2) mean(trans(B),2)]);
          disp(' ');
      end;
      if NF>0
          disp('Fixed coefficients');
          disp('RHOF    Accept rate since last printout');
          disp([RHOF kept./INFOSKIP]);
          disp('F');
          disp(F);
          disp(' ');
          kept=0;  %To reset after printout
      end
   end

   if r == NCREP
      disp('Burn-in over. Start saving draws.');
      disp(' ');
   end

   k=(r-NCREP)./NSKIP;
   if (r > NCREP) & (ceil(k) == (k))
       if NV>0
          adraws(:,k) = A;
          ddraws(:,:,k) = D;
          cmeans=cmeans+trans(B); %Divide by NEREP below to get mean instead of sum
       end
       if NF>0
          fdraws(:,k) = F;
       end
    end
end

if NV>0
   cmeans=cmeans./NEREP;
end
      

   