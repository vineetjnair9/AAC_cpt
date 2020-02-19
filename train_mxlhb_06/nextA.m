% Create next draw of A given D and B_n for all n
% Written by Kenneth Train on July 16, 2006


function anew=nextA

global B D NV SQRTNP IDV

anew=mean(B,2) + (chol(D)' * randn(NV,1))./SQRTNP;
anew(IDV(:,2) == 5)=0;