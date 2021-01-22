%% Running MSL with different hyperparams

clear; clc;
load('LL_finals.mat');
load('LL0s.mat');

% The distributions can be 
% 1. normal: N(b,w^2) where mean b and standard deviation w are estimated.
% 2. lognormal: coefficient is exp(beta) where beta~N(b,w^2) with b and w estimated
% 3. truncated normal, with the share below zero massed at zero: max(0,beta) where 
%                      beta~N(b,w^2) with b and w estimated.
% 4. S_B: exp(beta)/(1+exp(beta)) where beta~N(b,w^2) with b and w estimated.  [Johnson SB distribution?]
% 5. normal with zero mean (for error components): N(0,w^2) where w is estimated.
% 6. triangular: b+w*t where t is triangular between -1 and 1 and mean b and spread w are estimated.

% Type of draws to use in simulation
% 1=pseudo-random draws
% 2=standard Halton draws
% 3=shifted and shuffled Halton draws
% 4=modified Latin hypercube sampling, shifted and shuffled 
% 5=create your own draws or load draws from file

% SE = Inf occurs mainly due to no. of draws 
% But the critical no. of draws also depends on other factors
% Including draw type and starting values
run = 5;
draw_type = 4;
dist_type = 1;
b_init = -0.2;
w_init = 0.05;
init = 0.1;
ndraws = 500;
random = 1;

tic
if random
    mxlmsl_vineet_random;
else
    mxlmsl_vineet_fixed;
end
toc

%%
filename1 = join(["paramhat",num2str(run),".mat"],"");
save(filename1,'paramhat');

filename2 = join(["stderr",num2str(run),".mat"],"");
save(filename2,'stderr');

LL_finals(run) = LL_final;
LL0s(run) = LL0;
save('LL_finals.mat','LL_finals');
save('LL0s.mat','LL0s');