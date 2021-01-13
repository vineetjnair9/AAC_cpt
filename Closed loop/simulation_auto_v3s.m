%% Problem setup
% General reference case with gradient descent update (with sub-sampling ride offers)
3
% Consider a single SMoDS server operating in a particular locality/region
clear; clc;

% SMoDS server samples ride offers and decisions every dt seconds 
% Each sampling action also = one negotiation iteration
dt = 60;

% Total time horizon (s)
T = 3600; % (1h)

% Negotiation iterations (1 iteration = 1 new/updated/revised ride offer)
t = 1:dt:T; % Times at which we sample ride offers
num_iters = length(t);

% Time series of desired acceptance probabilities
% Mean = 0.5, SD = 0.1
p_star = (0.5 + 0.1.*randn(1,1)).*ones(1,num_iters);

% Update p_star after every delta iterations
delta = 360;
index = delta;
while index < num_iters
    p_star(index:end) = (0.5 + 0.1.*randn(1,1)).*ones(num_iters - index + 1,1);
    index = index + delta;
end

figure(1) 
plot(t,p_star)
xlabel('Time (s)','Interpreter','latex');
ylabel('Desired probability of acceptance $p^*(t)$','Interpreter','latex');
ylim([0 1]);

num_users = 500; % Total number of unique passengers using current server

%% True CPT parameters - Vector of parameters for each passenger
% 0 < alpha_plus, alpha_minus, beta_plus, beta_minus < 1, 1 < lambda < 100

theta_true = zeros(5,num_users);
theta_true(1:4,:) = rand(4,num_users);
theta_true(5,:) = 1 + 99.*rand(1,num_users);

%% Initializations for CPT params assumed by our model
theta_hat = zeros(5,num_users);
theta_hat(1:4,:) = rand(4,num_users);
theta_hat(5,:) = 1 + 99.*rand(1,num_users);

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
gamma_ub1 = 10;
gamma_ub2 = 30;

%% Construct trip characteristics
% Library of trips to use during iterations

num_trips = num_users; % Library size (# of trips/ride offer)

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

% Number of ride offer experiments to sample at each negotiation iteration
num_samples = 100; 

% Record actual p_sR values (averaged) to make convergence plot
p_sR_mean = zeros(1,num_iters);
p_sR_median = zeros(1,num_iters);

% Average & median absolute error at each iteration between p* and average p_sR value across all sampled ride offers
error_mean = zeros(1,num_iters);
error_median = zeros(1,num_iters);

% Interval for initial guess
x0_lb = (gamma_lb1 + gamma_lb2)*0.5;
x0_ub = (gamma_ub1 + gamma_ub2)*0.5;

% If we use fsolve
x0 = x0_lb + (x0_ub - x0_lb).*rand.*ones(1,num_samples);

options = optimoptions('fsolve','Display','off','MaxIterations',10000,'MaxFunctionEvaluations',100*num_samples,...
    'OptimalityTolerance',1e-4,'UseParallel',true);

% Initial SMoDS prices for each trip (random initialization)
gamma_curr = gamma_lb + (gamma_ub - gamma_lb).*rand(1,num_trips); % for current iteration/time step
gamma_prev = zeros(1,num_trips);
errors_curr = zeros(1,num_trips);
errors_prev = zeros(1,num_trips);
p_sR_true = zeros(1,num_trips);

R_type = 'dynamic_u0';
weight_type = 'original';
normalized = 1;
p_sR1 = @(gamma,u0S,x1S,x2S,pS,b_smS,alpha_plus,alpha_minus,beta_plus,beta_minus,lambda) accept_prob(gamma,u0S,x1S,x2S,pS,b_smS,R_type,weight_type,alpha_plus,alpha_minus,beta_plus,beta_minus,lambda);

tic
for i = 1:num_iters
            
    sampled_users = randsample(num_users,num_samples);
    
    % Vectors of parameters for sampled passengers
    alpha_plus_hat = theta_hat(1,sampled_users);
    alpha_minus_hat = theta_hat(2,sampled_users);
    beta_plus_hat = theta_hat(3,sampled_users);
    beta_minus_hat = theta_hat(4,sampled_users);
    lambda_hat = theta_hat(5,sampled_users);
    
    alpha_plus_true = theta_true(1,sampled_users);
    alpha_minus_true = theta_true(2,sampled_users);
    beta_plus_true = theta_true(3,sampled_users);
    beta_minus_true = theta_true(4,sampled_users);
    lambda_true = theta_true(5,sampled_users);
    
    % Sampled trip features
    x1S = x1(sampled_users);
    x2S = x2(sampled_users);
    pS = p(sampled_users);
    b_smS = b_sm(sampled_users);
    u0S = u0(sampled_users);
    
    p_sR2 = @(gamma,alpha_plus,alpha_minus,beta_plus,beta_minus,lambda) p_sR1(gamma,u0S,x1S,x2S,pS,b_smS,alpha_plus,alpha_minus,beta_plus,beta_minus,lambda);
    
    % Desired probability of acceptance at current iteration/sampling time (for all sampled users) 
    p_star_vecS = p_star(i).*ones(1,num_samples);
    p_star_vec = p_star(i).*ones(1,num_trips);
          
    % Ocptimal price calculation based on assumed CPT model parameters
    % Only do this to set the initial condition for each period/interval when p* changes
    if (i == 1) || (p_star(i) ~= p_star(i-1))

        % Numerical method: Solve f(gamma_hat) - p* = 0 
        
        % Compute subjective probability of acceptance
        p_sR_hat = @(gamma) p_sR2(gamma,alpha_plus_hat,alpha_minus_hat,beta_plus_hat,beta_minus_hat,lambda_hat);
        fun = @(gamma) p_sR_hat(gamma) - p_star_vecS;
        
        % Initialize solver at previous optimal tariff
        x0 = gamma_prev(sampled_users);
        
        % Solve nonlinear system of equations
        [gamma_currS,fval,exitflag,output] = fsolve(fun,x0,options);
               
        % Calculate actual probabilities of acceptance of passengers using their true parameters and the above price

        % Under assumption that true passenger follows assumed CPT model (same functional form) but with different parameters
        p_sR_trueS = p_sR2(gamma_currS,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        p_sR_true(sampled_users) = p_sR_trueS;
        
        % p_sR values calculated via CPT above will differ for each trip/ride offer being considered
        % Get single value by averaging across all the samples (for convergence plot)
        p_sR_mean(i) = mean(p_sR_true);
        p_sR_median(i) = median(p_sR_true);
        error_mean(i) = abs(p_sR_mean(i) - p_star(i));
        error_median(i) = abs(p_sR_median(i) - p_star(i));
        
        gamma_prev = gamma_curr;
        gamma_curr(sampled_users) = gamma_currS;
        
        errors_prev = errors_curr;
        errors_curr = p_star_vec - p_sR_true;
        
    % After initial condition, update tariff using only feedback until p_star changes
    else
        gamma_currS = gamma_curr(sampled_users); % Sampled subset
        gamma_prevS = gamma_prev(sampled_users);
        p_sR_trueS = p_sR2(gamma_currS,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        p_sR_true(sampled_users) = p_sR_trueS;
        p_sR_mean(i) = mean(p_sR_true);
        p_sR_median(i) = median(p_sR_true);
        error_mean(i) = abs(p_sR_mean(i) - p_star(i));
        error_median(i) = abs(p_sR_median(i) - p_star(i));
        
        errors_prev = errors_curr;
        errors_prevS = errors_prev(sampled_users);
        errors_curr = p_star_vec - p_sR_true;
        errors_currS = errors_curr(sampled_users);

        lr = 0.1; % Learning rate (step size) for gradient descent
        
        if normalized
            % Gradient descent update (output feedback signal to modify input)
            rel_delta_error = (errors_currS - errors_prevS)./errors_prevS;
            rel_delta_gamma = (gamma_currS - gamma_prevS)./gamma_prevS;

            % Normalized gradient
            gradS = rel_delta_error./rel_delta_gamma;

            % Set elements (trips) with gamma_currS = gamma_prevS to have zero graident (i.e. no update made)
            gradS(isnan(gradS)) = 0;
            
            % Using normalized gradient
            gamma_prev = gamma_curr;
            gamma_curr(sampled_users) = gamma_currS.*(1 - lr.*gradS);
        else
            % Gradient of error w.r.t tariff (non-normalized)
            gradS = (errors_currS - errors_prevS)./(gamma_currS - gamma_prevS);
            
            % Set elements (trips) with gamma_currS = gamma_prevS to have zero graident (i.e. no update made)
            gradS(isnan(gradS)) = 0;
            
            % Using non-normalized gradient
            gamma_prev = gamma_curr;
            gamma_curr(sampled_users) = gamma_currS - lr.*gradS;
        end
    end
         
end
toc

%% Plots
figure(1)
hold on
plot(t,p_star);
plot(t,p_sR_mean);
plot(t,p_sR_median);
ylim([0 1]);
xlabel('Time (s)','Interpreter','latex');
ylabel('Probability of acceptance','Interpreter','latex');
legend('Desired $p^*(t)$','Actual mean $p^s_R(t)$','Actual median $p^s_R(t)$','Interpreter','latex');

figure(2)
hold on
plot(t,error_mean);
plot(t,error_median);
legend('Mean','Median');
xlabel('Time (s)','Interpreter','latex');
ylabel('$|p^s_R - p^*|$','Interpreter','latex');
