%% Generating all possible combinations of MSL hyperparameters
% To use in 1 large for loop on supercloud

%% Case 1: All random coefficients (except 1 ASC)
% Create matrix of parameter combinations with columns: 
% (1) No. of draws
% (2) Draw type
% (3) Type of distribution
% (4) Starting values

clear;
clc;

ndraws_vals = 1000:200:2000;
draw_types = [1,2,3,4];
dist_types = [3,4]; % Dist. type for price and utility terms
b_inits = -0.2:0.1:0.2;
w_inits = 0.01:0.05:0.2;

sz = [length(ndraws_vals),length(draw_types),length(dist_types),length(b_inits),length(w_inits)];
n_runs = prod(sz);

params = zeros(n_runs,16); % All the b and w estimates
std_errors = zeros(n_runs,16);
LL_finals = zeros(n_runs,1);

for s = 1:n_runs
    [i,j,k,m,n] = ind2sub(sz,s);
    ndraws = ndraws_vals(i);
    draw_type = draw_types(j);
    dist_type = dist_types(k);
    b_init = b_inits(m);
    w_init= w_inits(m);
    
    % Using function method
    %     [params(s,:),std_errors(s,:),LL_final] = mxlmsl_vineet_random(ndraws,draw_type,dist_type,b_init,w_init);
    
    mxlmsl_vineet_random;
    params(s,:) = paramhat; % Need to negate price/time coefficients later on
    std_errors(s,:) = stderr;
    LL_finals(s) = LL_final;
end

%% Case 2: All fixed coefficients (except 1 ASC)
% Create matrix of parameter combinations with columns: 
% (1) No. of draws
% (2) Draw type
% (3) Starting values

clear;
clc;

ndraws_vals = 1000:200:2000;
draw_types = [1,2,3,4];
inits = -0.5:0.1:0;

sz = [length(ndraws_vals),length(draw_types),length(inits)];
n_runs = prod(sz);

params = zeros(n_runs,8); % All the coefficient estimates
std_errors = zeros(n_runs,8);
LL_finals = zeros(n_runs,1);

for s = 1:n_runs
    [i,j,k] = ind2sub(sz,s);
    ndraws = ndraws_vals(i);
    draw_type = draw_types(j);
    init = inits(k);
    mxlmsl_vineet_fixed;
    params(s,:) = paramhat; % Need to negate price/time coefficients later on
    std_errors(s,:) = stderr;
    LL_finals(s) = LL_final;
end


