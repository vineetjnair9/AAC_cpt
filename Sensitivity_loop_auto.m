%% Sensitivity analysis code
clear;
clc;

%% Nominal CPT parameter values 
%% Loop through several possible combinations of parameters
% param_combo = 0;
% for alpha = 0.04:0.05:0.99
%     for beta = 0.04:0.05:0.99
%         for lambda = 1:1:20
%             for p = 0.04:0.05:0.99
%                 param_combo = param_combo + 1; % unique tag to identify this sub-case 
%                 nominal = [alpha beta lambda p];
          
%% Loop through several possible scenarios
% Variables: u0, b_sm, ub, lb
u0_lb = -10;
u0_ub = 10;
b_sm_lb = -0.99;
b_sm_ub = -0.01;
p_lb1 = 1;
p_lb2 = 5;
p_ub1 = 5;
p_ub2 = 20;

num_scenarios = 100;
result.u0 = zeros(1,num_scenarios);
result.b_sm = zeros(1,num_scenarios);
result.x1 = zeros(1,num_scenarios);
result.x2 = zeros(1,num_scenarios);
result.valid = zeros(1,num_scenarios);

% tariff bounds
result.lb = zeros(1,num_scenarios); 
zeros(1,num_scenarios);

for i = 1:num_scenarios
    scenario = i;
    
    result.u0(i) = u0_lb + (u0_ub - u0_lb)*rand;
    result.b_sm(i) = b_sm_lb + (b_sm_ub - b_sm_lb)*rand;
   
    % Determine limiting values for x1 & x2 that ensure non-trivial choice scenario
    x1_ub = u0 - b_sm*lb;
    x2_lb = u0 - b_sm*ub;

    % Then randomly pick values
    result.x1(i) = x1_ub - abs(x1_ub)*rand; 
    result.x2(i) = x2_lb + abs(x2_lb)*rand;

    % Sanity check
    u1_max = result.x1(i) + result.b_sm(i)*lb;
    u2_min = result.x2(i) + result.b_sm(i)*ub;
    result.valid(i) = (u1_max <= u0 && u0 <= u2_min);

    scenario_vals.u0 = u0;
    scenario_vals.x1 = x1;
    scenario_vals.x2 = x2;
    scenario_vals.b_sm = b_sm;
    scenario_vals.R = 'Best case';
    scenario_vals.lb = lb;
    scenario_vals.ub = ub;

    theta = 'alpha';
    devs = linspace(-10,10,100); 

    [opt_soln_numeric,opt_cost_numeric,rel_sens_soln,rel_sens_cost,theta_vals,theta_nom] = derivations_numerical(theta,nominal,scenario_vals,devs,0);

