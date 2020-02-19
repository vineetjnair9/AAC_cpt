
%Take draws from the parameter space. 
%These draws are held in BETAS which is size NPxNVxNDRAWS
BETAS=randi(NGridPts,[NP,NV,NDRAWS]); %BETAS go from 1 to NGridPts
BETAS=(BETAS-1)./(NGridPts-1);   %Now BETAS go from zero to 1 inclusive
for r=1:NV
  BETAS(:,r,:)=COEF(r,1)+(COEF(r,2)-COEF(r,1)).*BETAS(:,r,:);  %Now BETAS go from lower limit to upper limit for each coefficient
end;