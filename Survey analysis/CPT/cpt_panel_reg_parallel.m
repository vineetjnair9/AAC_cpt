%% CPT parameters derivation
clear; clc;

% Optimize regularization hyperparameters
% Also draw the mode choice coefficients from specified distribution instead of just using mean values

% Mode-choice parameters
load('paramhat.mat');
load('valid_indices.mat');
T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_6-10_n955');
table = T.Fulllaunch610n955;

% Vary if table columns change
table = convertvars(table,[7:8,10:12,48:139],'double');
num_responses = length(valid_indices); % Total no. of survey responses

% Each row gives risk attitude parameters - alpha, beta+, beta-, lambda
cpt = zeros(num_responses,5);
error = zeros(num_responses,1); % Squared 2-norm of residual

weight_type = 'original';
cdf = true;
R_type = 'dynamic_SMoDS';

% Regularization parameters
reg = [sqrt(2); sqrt(2)]; % lambda1, lambda2

% Settings for non linear least squares
lb = [0,0,0,0,0];
ub = [1,1,1,1,1];
x0 = [0.5,0.5,0.5,0.5,0.5];

panel = @(x,respondent_num) panel_obj(x,respondent_num,table,R_type,weight_type,paramhat,cdf,reg,lb,ub);

options = optimoptions('lsqnonlin','Display','off','FiniteDifferenceStepSize',1e-2,'DiffMinChange',1e-2,'DiffMaxChange',0.1,...
    'MaxFunctionEvaluations',5000,'MaxIterations',1000);
% 'UseParallel',true,'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10,'FiniteDifferenceType','central'

%% Parallel (unable to run lsqnonlin in parallel, set 'UseParallel' = false)
tic
parfor i = 1:length(valid_indices)   
    j = valid_indices(i);
    fun = @(x) panel(x,j);
    [cpt(i,:),error(i)] = lsqnonlin(fun,x0,lb,ub,options);
end
toc

% R_type = dynamic_u0, dynamic_SMoDS, static_u0, static_SMoDS
% cdf = true (CPT), false (PT)
% weighting = original, modified
% optimization method = 'lsq' (nonlinear least squares), 'grid' (grid search), 'newton', 'global-opt' (Global optimization toolbox)

% With normalized lambda
cpt(:,5) = cpt(:,5) * 100;

