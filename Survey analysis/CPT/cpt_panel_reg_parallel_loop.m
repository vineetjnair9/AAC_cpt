clear; clc;

% Optimize regularization hyperparameters
% Also draw the mode choice coefficients from specified distribution instead of just using mean values

% Mode-choice parameters
load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Code\AAC_cpt\Survey analysis\Mode_choice\paramhat5.mat')
load('valid_indices.mat');

T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_6-10_n955');
table = T.Fulllaunch610n955;
table = convertvars(table,[7:8,10:12,48:139],'double');

weight_type = 'original';
cdf = true;
R_type = 'dynamic_SMoDS';

% Settings for non linear least squares
lb = [0,0,0,0,0.01];
ub = [1,1,1,1,1];
x0 = [0.5,0.5,0.5,0.5,0.5];

panel = @(x,respondent_num,reg) panel_obj(x,respondent_num,table,R_type,weight_type,paramhat,cdf,reg,lb,ub);

options = optimoptions('lsqnonlin','Display','off','FiniteDifferenceStepSize',1e-2,'DiffMinChange',1e-2,'DiffMaxChange',0.1,...
    'MaxFunctionEvaluations',10000,'MaxIterations',10000,'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'FiniteDifferenceType','central');
% 'UseParallel',true

num_responses = length(valid_indices); % Total no. of survey responses
% Regularization parameters
reg_vals = sqrt(0:0.5:5);
sz = [num_responses,length(reg_vals)];
n_runs = prod(sz);

% Each row gives: [alpha+, alpha-, beta+, beta-, lambda, error]
% cpt = zeros(n_runs,5);
% error = zeros(n_runs,1);
% info = zeros(n_runs,2);
% row = [respondent, reg^2 = \nu, cpt_params(5), squared error_norm]
main = zeros(n_runs,8);

%% Parallel (unable to run lsqnonlin in parallel, set 'UseParallel' = false)
tic
for i = 1:n_runs  
    [j,k] = ind2sub(sz,i);
    respondent = valid_indices(j);
    reg_val = reg_vals(k);
    fun = @(x) panel(x,respondent,reg_val);
    main(i,1) = respondent;
    main(i,2) = reg_val^2;
    [main(i,3:7),main(i,8)] = lsqnonlin(fun,x0,lb,ub,options);
end
toc

% R_type = dynamic_u0, dynamic_SMoDS, static_u0, static_SMoDS
% cdf = true (CPT), false (PT)
% weighting = original, modified
% optimization method = 'lsq' (nonlinear least squares), 'grid' (grid search), 'newton', 'global-opt' (Global optimization toolbox)

% With normalized lambda
main(:,7) = main(:,7) * 100;
run = 'SMoDS_reg_OGweight_cdf_R-dynSMoDS_norm-trueCE_local.mat';
save(run,'main');
