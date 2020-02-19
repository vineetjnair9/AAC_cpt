% Creates next draw of D given A and B when D is constrained to be diagonal
% (ie, no correlations among coefficients.)
% Written by Kenneth Train on July 16, 2006.

function dnew=nextDdiag

global A B NV NP

s=B-repmat(A,1,NP);
dnew=1+sum(s.^2,2);
s=repmat(1./sqrt(dnew),1,NP+1).*randn(NV,NP+1);
dnew=1./sum(s.^2,2);
dnew=diag(dnew);