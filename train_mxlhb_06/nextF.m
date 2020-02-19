% Create next draw of F 
% Written by Kenneth Train on July 25, 2006.

function [fnew,pnew,ind]=nextF;

global B F P RHOF NF


fnew= F + sqrt(RHOF) .* randn(NF,1);


pnew=logit(fnew,trans(B));

r=prod((pnew ./ P),2);
ind= (rand(1,1) <= r);


pnew=pnew.*ind + P.*(1-ind);
fnew=fnew.*ind + F.*(1-ind);
