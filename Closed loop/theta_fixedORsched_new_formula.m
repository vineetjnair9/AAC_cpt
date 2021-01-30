%% Problem setup
% General reference case with gradient descent update (without sub-sampling ride offers)
% Fixed step-size and learning rate schedules

% Consider a single SMoDS server operating in a particular locality/region
clear; clc;

lr_fixed = 0.5; 
time_factor = 0.005;
exp_factor = 0.005;
step_factor = 0.005;
step_iters = 5;
lr_type = 'fixed';  % 'fixed', 'time-decay', 'exp-decay', 'step-decay'
R_type = 'dynamic_SMoDS';
weight_type = 'original';

% SMoDS server samples ride offers and decisions every dt seconds 
% Each sampling action also = one negotiation iteration
dt = 30;

% Total time horizon (s)
T = 1800; % (1h)

% Negotiation iterations (1 iteration = 1 new/updated/revised ride offer)
t = 1:dt:T; % Times at which we sample ride offers
num_iters = length(t);

% Time series of desired acceptance probabilities
Mean = 0.5; SD = 0.1;
p_star = (Mean + SD.*randn(1,1)).*ones(1,num_iters);

% Update p_star after every delta iterations
delta = round(num_iters/1);
index = delta;
while index < num_iters
    p_star(index:end) = (0.5 + 0.1.*randn(1,1)).*ones(num_iters - index + 1,1);
    index = index + delta;
end

num_users = 500; % Total number of unique passengers using current server

% True CPT parameters - Vector of parameters for each passenger
% 0 < alpha_plus, alpha_minus, beta_plus, beta_minus < 1, 1 < lambda < 100

theta_true = zeros(5,num_users);
theta_true(1:4,:) = rand(4,num_users);
theta_true(5,:) = 1 + 99.*rand(1,num_users);

% Initializations for CPT params assumed by our model
theta_hat = zeros(5,num_users);
theta_hat(1:4,:) = rand(4,num_users);
theta_hat(5,:) = 1 + 99.*rand(1,num_users);

% Bounds on range from which to draw parameters uniformly at random

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
gamma_ub1 = 15;
gamma_ub2 = 30;

% Construct trip characteristics
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

% Ride offer iterations
% Record actual p_sR values (averaged) to make convergence plot
p_sR_mean = zeros(1,num_iters);

% Squared norm of error (objective function)
error_norm = zeros(1,num_iters);
error_mean = zeros(1,num_iters);

% Interval for initial guess
x0_lb = (gamma_lb1 + gamma_lb2)*0.5;
x0_ub = (gamma_ub1 + gamma_ub2)*0.5;

% If we use fsolve
x0 = x0_lb + (x0_ub - x0_lb).*rand.*ones(1,num_trips);

options = optimoptions('fsolve','Display','off','MaxIterations',10000,'MaxFunctionEvaluations',100*num_trips,'UseParallel',true);

% Initial SMoDS prices for each trip (random initialization)
gamma_curr = gamma_lb + (gamma_ub - gamma_lb).*rand(1,num_trips); % for current iteration/time step
p_sR_curr = zeros(1,num_trips);
p_sR_prev = zeros(1,num_trips);

% Vectors of parameters
alpha_plus_curr = theta_hat(1,:);
alpha_minus_curr = theta_hat(2,:);
beta_plus_curr = theta_hat(3,:);
beta_minus_curr = theta_hat(4,:);
lambda_curr = theta_hat(5,:);

alpha_plus_true = theta_true(1,:);
alpha_minus_true = theta_true(2,:);
beta_plus_true = theta_true(3,:);
beta_minus_true = theta_true(4,:);
lambda_true = theta_true(5,:);

p_sR1 = @(gamma,alpha_plus,alpha_minus,beta_plus,beta_minus,lambda) accept_prob(gamma,u0,x1,x2,p,b_sm,R_type,weight_type,alpha_plus,alpha_minus,beta_plus,beta_minus,lambda);

tic
for i = 1:num_iters

    % Desired probability of acceptance at current iteration/sampling time
    p_star_vec = p_star(i).*ones(1,num_trips);
          
    % Optimal price calculation based on assumed CPT model parameters

    if (i == 1) % Need to randomly perturb the CPT parameters to get data for 1st iter
        alpha_plus_prev = alpha_plus_curr;
        alpha_minus_prev = alpha_minus_curr;
        beta_plus_prev = beta_plus_curr;
        beta_minus_prev = beta_minus_curr;
        lambda_prev = lambda_curr;
        
        alpha_plus_curr = alpha_plus_prev.*rand(1,num_trips);
        alpha_minus_curr = alpha_minus_prev.*rand(1,num_trips);
        beta_plus_curr = beta_plus_prev.*rand(1,num_trips);
        beta_minus_curr = beta_minus_prev.*rand(1,num_trips);
        lambda_curr = lambda_prev.*rand(1,num_trips);
        
        % Numerical method: Solve f(gamma_hat) - p* = 0 
        % Compute subjective probability of acceptance
        p_sR_hat = @(gamma) p_sR1(gamma,alpha_plus_curr,alpha_minus_curr,beta_plus_curr,beta_minus_curr,lambda_curr);
        fun = @(gamma) p_sR_hat(gamma) - p_star_vec;
                
        % Solve nonlinear system of equations
        gamma_prev = gamma_curr;
        gamma_curr = fsolve(fun,x0,options);
               
        % Calculate actual probabilities of acceptance of passengers using their true parameters and the above price
        % Under assumption that true passenger follows assumed CPT model (same functional form) but with different parameters
        p_sR_prev = p_sR_curr;
        p_sR_curr = p_sR1(gamma_curr,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        
        % p_sR values calculated via CPT above will differ for each trip/ride offer being considered
        % Get single value by averaging across all the samples (for convergence plot)
        p_sR_mean(i) = mean(p_sR_curr);
        error_norm(i) = 0.5*norm(p_sR_curr - p_star(i))^2;
        error_mean(i) = abs(p_sR_mean(i) - p_star(i));

    else    
        x0 = gamma_prev;
        % Gradient descent update (output feedback signal to modify input)
        % Learning rate (step size) for gradient descent
        if strcmp(lr_type,'fixed')
            % Fixed learning rate 
            lr = lr_fixed;
        elseif strcmp(lr_type,'time-decay')
            % Learning rate schedules 
            % Time-based decay
            lr = lr_fixed/(1 + time_factor*n);
        elseif strcmp(lr_type,'exp-decay')
            % Exponential decay
            lr = lr_fixed*exp(-exp_factor*n);
        else
            lr = lr_fixed*step_factor^(floor(n/step_iters));
        end
            
        % Gradient of error w.r.t tariff (non-normalized)
        grad_alpha_plus = (p_sR_curr - p_sR_prev)./(alpha_plus_curr - alpha_plus_prev);
        grad_alpha_minus = (p_sR_curr - p_sR_prev)./(alpha_minus_curr - alpha_minus_prev);
        grad_beta_plus = (p_sR_curr - p_sR_prev)./(beta_plus_curr - beta_plus_prev);
        grad_beta_minus = (p_sR_curr - p_sR_prev)./(beta_minus_curr - beta_minus_prev);
        grad_lambda = (p_sR_curr - p_sR_prev)./(beta_minus_curr - beta_minus_prev);
        
        % Set elements (trips) with gamma_currS = gamma_prevS to have zero graident (i.e. no update made)
        grad_alpha_plus(isnan(grad_alpha_plus)) = 0;
        grad_alpha_minus(isnan(grad_alpha_minus)) = 0;
        grad_beta_plus(isnan(grad_beta_plus)) = 0;
        grad_beta_minus(isnan(grad_beta_minus)) = 0;
        grad_lambda(isnan(grad_lambda)) = 0;
        
        update_alpha_plus = (p_sR_curr - p_star_vec).*grad_alpha_plus;
        update_alpha_minus = (p_sR_curr - p_star_vec).*grad_alpha_minus;
        update_beta_plus = (p_sR_curr - p_star_vec).*grad_beta_plus;
        update_beta_minus = (p_sR_curr - p_star_vec).*grad_beta_minus;
        update_lambda = (p_sR_curr - p_star_vec).*grad_lambda;
        
        alpha_plus_prev = alpha_plus_curr;
        alpha_minus_prev = alpha_minus_curr;
        beta_plus_prev = beta_plus_curr;
        beta_minus_prev = beta_minus_curr;
        lambda_prev = lambda_curr;
        
        alpha_plus_curr = alpha_plus_curr - lr.*update_alpha_plus;
        alpha_minus_curr = alpha_minus_curr - lr.*update_alpha_minus;
        beta_plus_curr = beta_plus_curr	- lr.*update_beta_plus;
        beta_minus_curr = beta_minus_curr - lr.*update_beta_minus;
        lambda_curr = lambda_curr - lr.*update_lambda;
 
        p_sR_hat = @(gamma) p_sR1(gamma,alpha_plus_curr,alpha_minus_curr,beta_plus_curr,beta_minus_curr,lambda_curr);
        fun = @(gamma) p_sR_hat(gamma) - p_star_vec;
                
        % Solve nonlinear system of equations
        gamma_prev = gamma_curr;
        gamma_curr = fsolve(fun,x0,options);
               
        % Calculate actual probabilities of acceptance of passengers using their true parameters and the above price
        % Under assumption that true passenger follows assumed CPT model (same functional form) but with different parameters
        p_sR_prev = p_sR_curr;
        p_sR_curr = p_sR1(gamma_curr,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        
        % Get single value by averaging across all the samples (for convergence plot)
        p_sR_mean(i) = mean(p_sR_curr);
        error_norm(i) = 0.5*norm(p_sR_curr - p_star(i))^2;
        error_mean(i) = abs(p_sR_mean(i) - p_star(i));
    end
         
end
toc

% Plots
figure(1)
hold on
plot(t,p_star,'LineWidth',1.2);
plot(t,p_sR_mean,'LineWidth',1.2);
ylim([0 1]);
xlabel('Time (s)','Interpreter','latex');
ylabel('Probability of acceptance','Interpreter','latex');
legend('Desired $p^*(t)$','Actual mean $p^s_R(t)$','Interpreter','latex');

figure(2)
yyaxis left
plot(t,error_norm,'LineWidth',1.2);
xlabel('Time (s)','Interpreter','latex');
ylabel('$\frac{1}{2}\|p^s_R - p^*\|^2$','Interpreter','latex');
yyaxis right
plot(t,error_mean,'LineWidth',1.2);
ylabel('Mean $|p^s_R - p^*|$','Interpreter','latex');
