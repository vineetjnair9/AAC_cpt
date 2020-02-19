% Creates next draw of D given A and B.
% Full covariance D. Use nextDdiag for diagnonal D (ie, uncorrelated coefficients.)
% Written by Kenneth Train on July 16, 2006.

function dnew=nextD

global A B NV NP
 
s=B-repmat(A,1,NP);
dnew= (s*s') + NV.*eye(NV);
s=chol(inv(dnew))' * randn(NV,NP+NV);
dnew=inv(s*s');

