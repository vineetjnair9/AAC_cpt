% Transform normally distributed terms into coefficients
% Written by Kenneth Train, July 14, 2006.

% Input b has dimension NVxNP. 
% Output c has dimension NVxNP.
% Uses IDV to determine transformations.


function c=trans(b)
 
 global IDV NV

 c=b;
 if NV>0
    c(IDV(:,2) == 2,:)=exp(c(IDV(:,2) == 2,:));
    c(IDV(:,2) == 3,:)=c(IDV(:,2) == 3,:).*(c(IDV(:,2) == 3,:)>0);
    c(IDV(:,2) == 4,:)=exp(c(IDV(:,2) == 4,:))./(1+exp(c(IDV(:,2) == 4,:)));
 end