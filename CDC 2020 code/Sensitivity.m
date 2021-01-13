%% Sensitivity analysis code
clear;
clc;

%% Nominal CPT parameter values 
% Taken from CPT CDC paper
alpha = 0.895;
beta = 0.88;
lambda = 2.25;
p = 0.7;
nominal = [alpha beta lambda p];

%% Construct scenarios
% Method 1: Generate scenarios by actually considering terms in utility
% function individually
a_sm = [-0.2 -0.1 -0.09]';
b_sm = -0.18;
c_sm = -1.7;

% Using any other alternative (certain prospect)
a_0 = [-0.06 -0.04 -0.08]';
b_0 = -0.14;
c_0 = -1.3;
gamma_0 = 14.79; % ($)
ub = 5;
lb = 2.5;

% Travel times (in min)
t0 = [0 4 10]'; 
t1 = [5.81 10.96 19.21]';
t2 = [5.81 4.96 15.21]'; 

u0 = a_0'*t0 + b_0*gamma_0 + c_0;
x1 = a_sm'*t1;
x2 = a_sm'*t2;

scenario_vals.b_sm = b_sm;
scenario_vals.b_0 = b_0;
scenario_vals.R = 'Best case';
scenario_vals.lb = lb;
scenario_vals.ub = ub;
scenario_vals.u0 = u0;
scenario_vals.x1 = x1;
scenario_vals.x2 = x2;

theta = 'alpha';
devs = linspace(-20,20,100); % Consider parameter perturbations between +/- 5%

[dGamma_dTheta_nom,dF_dTheta_nom,error,max_dev] = derivations_v2(theta,nominal,scenario_vals,devs);

%% Method 2: Just vary u0, x1 and x2

u0 = -4;
b_sm = -0.18;
ub = 5;
lb = 2.5;

% Determine limiting values for x1 & x2 that ensure non-trivial choice scenario
x1_ub = u0 - b_sm*lb;
x2_lb = u0 - b_sm*ub;

% Then randomly pick values within +/- 20%
x1 = x1_ub - abs(x1_ub)*0.2*rand;
x2 = x2_lb + abs(x2_lb)*0.2*rand;

% Sanity check
u1_max = x1 + b_sm*lb;
u2_min = x2 + b_sm*ub;
Check = (u1_max <= u0 && u0 <= u2_min);

scenario_vals.u0 = u0;
scenario_vals.x1 = x1;
scenario_vals.x2 = x2;
scenario_vals.b_sm = b_sm;
scenario_vals.R = 'Best case';
scenario_vals.lb = lb;
scenario_vals.ub = ub;

theta = 'alpha';
devs = linspace(-20,20,100); % Consider parameter perturbations between +/- 5%

[dGamma_dTheta_nom,dF_dTheta_nom,opt_soln_analytic,opt_cost_analytic,error,max_dev,rel_sens_soln,rel_sens_cost_actual] = derivations_local(theta,nominal,scenario_vals,devs,0);
[opt_soln_numeric,opt_cost_numeric,rel_sens_soln,rel_sens_cost,mismatch_loss,theta_vals,theta_nom] = derivations_numerical(theta,nominal,scenario_vals,devs,0);

%% Comparison between analytical and numerical results

% Compute relative error treating numerical solution as true value

% Error in 1st order approximation of optimal solution
Opt_soln_error = (norm(opt_soln_analytic - opt_soln_numeric)/norm(opt_soln_numeric))*100

Opt_cost_error = (norm(opt_cost_analytic + opt_cost_numeric)/norm(opt_cost_numeric))*100

figure(1)

subplot(2,1,1)
hold on 
plot(theta_vals,opt_soln_analytic);
plot(theta_vals,opt_soln_numeric);
xlabel('Parameter value');
ylabel('Optimal value of objective function ($)')    
legend('Analytical','Numerical');
hold off

subplot(2,1,2)
hold on 
plot(theta_vals,opt_cost_analytic);
plot(theta_vals,-opt_cost_numeric);
xlabel('Parameter value');
ylabel('Optimal value of objective function ($)')    
legend('Analytical','Numerical');
hold off
