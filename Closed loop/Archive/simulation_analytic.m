%% Problem setup
% Using analytically derived acceptance probability formula

% Consider a single SMoDS server operating in a particular locality/region
clear; clc;

% Set random seed for reproduceability of results
% rng('default')

% SMoDS server samples ride offers and decisions every dt seconds 
% Each sampling action also = one negotiation iteration
dt = 10;

% Number of ride offer experiments to sample at each negotiation iteration
num_samples = 500; 

% Total time horizon (s)
T = 3600; % (1h)

% Negotiation iterations
t = 1:dt:T; % Times at which
num_iters = length(t);

% Time series of desired acceptance probabilities
mu = 0.5; sd = 0.15; 
p_star = (mu + sd.*randn(1,1)).*ones(1,num_iters);

% Update p_star after every delta iterations
delta = 100;
index = delta;
while index < num_iters
    p_star(index:end) = (mu + sd.*randn(1,1)).*ones(num_iters - index + 1,1);
    index = index + delta;
end

figure(1) 
plot(t,p_star)
xlabel('Time (s)','Interpreter','latex');
ylabel('Desired probability of acceptance $p^*(t)$','Interpreter','latex');
ylim([0 1]);

num_users = 500; % Total number of unique passengers using current server

%% True CPT parameters - Vector of parameters for each passenger
% To simplify, assume alpha_plus = alpha_minus, beta_plus = beta_minus
% 1 < alpha, beta < 1, 1 < lambda < 100

theta_true = zeros(3,num_users);
theta_true(1:2,:) = rand(2,num_users);
theta_true(3,:) = 1 + 99.*rand(1,num_users);

%% Initializations for CPT params assumed by our model
theta_hat = zeros(3,num_users);
theta_hat(1:2,:) = rand(2,num_users);
theta_hat(3,:) = 1 + 99.*rand(1,num_users);

%% Bounds on range from which to draw parameters uniformly at random

% Objective utility of certain travel alternative 
% e.g. exclusive ridehailing or public transit
u0_lb = -10;
u0_ub = 10;

% Disutility due to price for SMoDS
b_sm_lb = -0.99;
b_sm_ub = -0.01;

% For this problem, assume both SMoDS outcomes have same price (differ only in travel times)
% In order to make price-based control easier/more tractable
% Lower and upper bounds on price for outcomes u1 and u2
gamma_lb1 = 3;
gamma_lb2 = 5;
gamma_ub1 = 8;
gamma_ub2 = 20;

%% Construct trip characteristics
% Library of trips to use during iterations

num_trips = num_users; % Library size (# of trips)

u0 = u0_lb + (u0_ub - u0_lb).*rand(1,num_trips);
b_sm = b_sm_lb + (b_sm_ub - b_sm_lb).*rand(1,num_trips);
gamma_lb = gamma_lb1 + (gamma_lb2 - gamma_lb1).*rand(1,num_trips);
gamma_ub = gamma_ub1 + (gamma_ub2 - gamma_ub1).*rand(1,num_trips);

% u1 --> u_low, u2 --> u_high
% Determine limiting values for x1 & x2 that ensure non-trivial choice scenario
x1_ub = u0 - b_sm.*gamma_lb;
x2_lb = u0 - b_sm.*gamma_ub;

% Then randomly pick values
x1 = x1_ub - abs(0.8.*x1_ub).*rand; 
x2 = x2_lb + abs(0.8.*x2_lb).*rand;

% Sanity check
u1_max = x1 + b_sm.*gamma_lb;
u2_min = x2 + b_sm.*gamma_ub;
valid = sum(u1_max <= u0 & u0 <= u2_min) == num_trips;

% Model SMoDS as Bernoulli dist with only 2 outcomes (u1,u2) ~ (p,1-p)
p = rand(1,num_trips);

%% Ride offer iterations
% Case 1: Using gamma directly as the control input being varied

% Record actual p_sR values (averaged) to make convergence plot
p_sR_vals = zeros(1,num_iters);

% Average absolute error at each iteration between p* and average p_sR value across all sampled ride offers
error_vals = zeros(1,num_iters);

% Interval for initial guess
x0_lb = (gamma_lb1 + gamma_lb2)*0.5;
x0_ub = (gamma_ub1 + gamma_ub2)*0.5;

% % If we use fzero
% x0 = [x0_lb x0_ub];

% If we use fsolve
x0 = x0_lb + (x0_ub - x0_lb).*rand.*ones(1,num_samples);

options = optimoptions('fsolve','Display','off','MaxIterations',10000,'MaxFunctionEvaluations',100*num_samples,...
    'OptimalityTolerance',1e-4,'UseParallel',true);

% Initial SMoDS prices for each trip
gamma_prev = gamma_lb + (gamma_ub - gamma_lb).*rand(1,num_trips);

tic
for i = 1:num_iters
        
    sampled_users = randsample(num_users,num_samples);
    
    % Vectors of parameters for sampled passengers
    alpha_hat = theta_hat(1,sampled_users);
    beta_hat = theta_hat(2,sampled_users);
    lambda_hat = theta_hat(3,sampled_users);
    alpha_true = theta_true(1,sampled_users);
    beta_true = theta_true(2,sampled_users);
    lambda_true = theta_true(3,sampled_users);
    
    % Sampled trip features
    x1S = x1(sampled_users);
    x2S = x2(sampled_users);
    pS = p(sampled_users);
    b_smS = b_sm(sampled_users);
    u0S = u0(sampled_users);
    
    % Desired probability of acceptance at current iteration/sampling time (for all sampled users) 
    p_star_vec = p_star(i).*ones(1,num_samples);
          
    % Ocptimal price calculation based on assumed CPT model parameters
    % Only do this to set the initial condition for each period/interval when p* changes
    if (i == 1) || (p_star(i) ~= p_star(i-1))
    %     % Analytical method: Using symbolic calculations 
    %     syms Gamma [1 num_samples]
    %     p_sR = @(Gamma) accept_prob(Gamma,alpha_hat,beta_hat,lambda_hat,x1,x2,p,b_sm,u0);
    %     gamma_hat = finverse(p_sR); 
    %     gamma_hat_vals = subs(gamma_hat,p_stars);

        % Numerical method: Solve f(gamma_hat) - p* = 0 
        p_sR = @(gamma) accept_prob_analytic(gamma,alpha_hat,beta_hat,lambda_hat,x1S,x2S,pS,b_smS,u0S);
        fun = @(gamma) p_sR(gamma) - p_star_vec;
        
        % Initialize solver at previous optimal tariff
        x0 = gamma_prev(sampled_users);
        
        % Solve nonlinear system of equations
        [Gamma,fval,exitflag,output] = fsolve(fun,x0,options);
        
%         % With fzero, need to loop through each point since it can only be used for scalar valued functions
%         Gamma = fzero(fun,x0); 
        
        % Calculate actual probabilities of acceptance of passengers using their true parameters and the above price

        % Under assumption that true passenger follows assumed CPT model (same functional form) but with different parameters
        p_sR = accept_prob_analytic(Gamma,alpha_true,beta_true,lambda_true,x1S,x2S,pS,b_smS,u0S);
        
        % p_sR values calculated via CPT above will differ for each trip/ride offer being considered
        % Get single value by averaging across all the samples (for convergence plot)
        p_sR_vals(i) = mean(p_sR);
        error_vals(i) = abs(p_sR_vals(i) - p_star(i));

        % Note: Could also use empirical probabilities (calculated numerically) if we have access to actual data
        % Observed decisions by passengers at this iteration (negotiation round)
        % Number of successes (passenger acceptances) in num_samples trials (ride offers)
        % Draw from binomial distribution: B(num_samples,p*)
        
    % After initial condition, update tariff using only feedback until p_star changes
    else
%         % Assume true reference function/mapping is known
%         % For now, set R = best case SMoDS outcome (same as CDC paper)
%         R = u2;

        gamma_prevS = gamma_prev(sampled_users);
        p_sR = accept_prob_analytic(gamma_prevS,alpha_true,beta_true,lambda_true,x1S,x2S,pS,b_smS,u0S);
        p_sR_vals(i) = mean(p_sR);
        error_vals(i) = abs(p_sR_vals(i) - p_star(i));
        
        errors =  p_star_vec - p_sR;

        % Control law
        K = 10;
        Gamma = gamma_prevS - K.*errors;  
    end
       
    gamma_prev(sampled_users) = Gamma;
end
toc

%% Plots
figure(1)
hold on
plot(t,p_star);
plot(t,p_sR_vals);
ylim([0 1]);
xlabel('Time (s)','Interpreter','latex');
ylabel('Probability of acceptance','Interpreter','latex');
legend('Desired $p^*(t)$','Actual $p^s_R(t)$','Interpreter','latex');

figure(2)
plot(t,error_vals);
xlabel('Time (s)','Interpreter','latex');
ylabel('$|p^s_R - p^*|$','Interpreter','latex');

% Also plot prices here for a few selected users/trips, varying over time/iterations
