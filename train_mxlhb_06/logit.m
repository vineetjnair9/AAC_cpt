% Calculate logit probability for chosen alternatives for each person
%    with multiple choice situations for each person and multiple people,
%    using globals for all inputs except coefficients.
% Written by Kenneth Train, first version on July 14, 2006..
%    lastest revision on April 3, 2007.
%
% Probability is Prod_t { exp(V_*t)/SUM_j[exp(V_jt)] } 
%             = Prod_t { 1 / (1+ Sum_j~=* [exp(V_jt-V-*t)] }
%    where * denotes chosen alternative, j is alternative, t is choice situation.
% Using differences from chosen alternative reduces computation time (with
%    one less alternative), eliminates need to retain and use the dependent
%    variable, and avoids numerical problems when exp(V) is very large or
%    small, since prob is 1/(1+k) which can be evaluated for any k, including
%    infinite k and infinitesimally small k. In contrast. e(V)/sum e(V) can
%    result in numerical "divide by zero" if denominator is sufficiently small
%    and NaN if numerator and denominator are both numerically zero.
% 
% Input f contains the fixed coefficients, and has dimension NFX1.
% Input c contains the random coefficients for each person, and has dimensions NV x NP.
% Either input can be an empty matrix. 
% Output p contains the logit probabilities, which is a row vector of dimension 1xNP.
% Code assumes that all GLOBALS already exist.


function p=logit(f,c);

global NV NP NALTMAX NCSMAX X S XF NF


if NF > 0
   ff=reshape(f,1,1,NF,1);
   ff=repmat(ff,[NALTMAX-1,NCSMAX,1,NP]);
   v=reshape(sum(XF.*ff,3),NALTMAX-1,NCSMAX,NP);  %w is (NALTMAX-1) x NCSMAX x NP
else
   v=zeros(NALTMAX-1,NCSMAX,NP);
end

if NV >0
   cc=reshape(c,1,1,NV,NP);
   cc=repmat(cc,[NALTMAX-1,NCSMAX,1,1]);
   v=v+reshape(sum(X.*cc,3),NALTMAX-1,NCSMAX,NP); %v is (NALTMAX-1)xNCSMAXxNP
end

v=exp(v);
v=1./(1+sum(v.*S,1)); %Logit prob for each choice situation 1xNCSMAXxNP
v=reshape(v,NCSMAX,NP);  %v is now NCSMAXxNP
w=reshape(sum(S,1),NCSMAX,NP);
v(w == 0)=1; %Make prob=1 for choice situations that person did not face any alts
p=prod(v,1); %p is 1xNP
p(1,isnan(p))=0; %Change missing values to 0, as a precaution.