% Create next draw of B_n for all n, given A and D
% Written by Kenneth Train on July 16, 2006.

function [bnew,pnew,newrho,i]=nextB;

global A B D P RHO NP NV F


bnew= B + sqrt(RHO) .* (chol(D)'*randn(NV,NP));

bn=bnew-repmat(A,1,NP);
bo=B-repmat(A,1,NP);

pnew=logit(F,trans(bnew)); 

r=(pnew ./ P) .* exp(-0.5 .* (sum(bn.*(inv(D)*bn)) - sum(bo.*(inv(D)*bo))));
ind= (rand(1,NP) <= r);
i=sum(ind')./NP;
newrho=RHO - 0.001.*(i < 0.3) + 0.001.*(i > 0.3);

pnew=pnew.*ind + P.*(1-ind);
bnew=bnew.*repmat(ind,NV,1) + B.*repmat((1-ind),NV,1);
