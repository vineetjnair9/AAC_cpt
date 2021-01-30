%% Problem setup
% General reference case with gradient descent update (without sub-sampling ride offers)
% Fixed step-size and learning rate schedules

% Consider a single SMoDS server operating in a particular locality/region
clear; clc;

lr_fixed = 0.3; 
time_factor = 0.005;
exp_factor = 0.005;
step_factor = 0.005;
step_iters = 5;
R_type = 'dynamic_SMoDS';
weight_type = 'original';
num_users = 1; % Total number of unique passengers using current server

% SMoDS server samples ride offers and decisions every dt seconds 
% Each sampling action also = one negotiation iteration
dt = 10;

% Total time horizon (s)
T = 1800; % 

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
gamma_ub1 = 10;
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
% If using purely random p
p = rand(1,num_trips);

% Ride offer iterations
% Case 1: Using gamma directly as the control input being varied

% Record actual p_sR values (averaged) to make convergence plot
p_sR_mean_fixed = zeros(1,num_iters);
p_sR_mean_time = zeros(1,num_iters);
p_sR_mean_exp = zeros(1,num_iters);
p_sR_mean_step = zeros(1,num_iters);

% Average & median absolute error at each iteration between p* and average p_sR value across all sampled ride offers
error_mean_fixed = zeros(1,num_iters);
error_mean_time = zeros(1,num_iters);
error_mean_exp = zeros(1,num_iters);
error_mean_step = zeros(1,num_iters);

error_norm_fixed = zeros(1,num_iters);
error_norm_time = zeros(1,num_iters);
error_norm_exp = zeros(1,num_iters);
error_norm_step = zeros(1,num_iters);

% Interval for initial guess
x0_lb = (gamma_lb1 + gamma_lb2)*0.5;
x0_ub = (gamma_ub1 + gamma_ub2)*0.5;

% If we use fsolve
x0_fixed = x0_lb + (x0_ub - x0_lb).*rand.*ones(1,num_trips);
x0_time = x0_fixed;
x0_exp = x0_fixed;
x0_step = x0_fixed;

options = optimoptions('fsolve','Display','off','MaxIterations',10000,'MaxFunctionEvaluations',100*num_trips,'UseParallel',true);

% Initial SMoDS prices for each trip (random initialization)
gamma_curr_fixed = gamma_lb + (gamma_ub - gamma_lb).*rand(1,num_trips); % for current iteration/time step
gamma_curr_time = gamma_curr_fixed;
gamma_curr_exp = gamma_curr_fixed;
gamma_curr_step = gamma_curr_fixed;

gamma_prev_fixed = zeros(1,num_trips);
gamma_prev_time = zeros(1,num_trips);
gamma_prev_exp = zeros(1,num_trips);
gamma_prev_step = zeros(1,num_trips);

p_sR_curr_fixed = zeros(1,num_trips);
p_sR_curr_time = zeros(1,num_trips);
p_sR_curr_exp = zeros(1,num_trips);
p_sR_curr_step = zeros(1,num_trips);

p_sR_prev_fixed = zeros(1,num_trips);
p_sR_prev_time = zeros(1,num_trips);
p_sR_prev_exp = zeros(1,num_trips);
p_sR_prev_step = zeros(1,num_trips);

% Vectors of parameters
alpha_plus_hat = theta_hat(1,:);
alpha_minus_hat = theta_hat(2,:);
beta_plus_hat = theta_hat(3,:);
beta_minus_hat = theta_hat(4,:);
lambda_hat = theta_hat(5,:);

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
          
    % Ocptimal price calculation based on assumed CPT model parameters
    % Only do this to set the initial condition for each period/interval when p* changes
    if (i == 1) || (p_star(i) ~= p_star(i-1))
        n = 0;

        % Numerical method: Solve f(gamma_hat) - p* = 0 
        % Compute subjective probability of acceptance
        p_sR_hat = @(gamma) p_sR1(gamma,alpha_plus_hat,alpha_minus_hat,beta_plus_hat,beta_minus_hat,lambda_hat);
        fun = @(gamma) p_sR_hat(gamma) - p_star_vec;
        
        % Initialize solver at previous optimal tariff
        if (i ~= 1)
            x0_fixed = gamma_prev_fixed;
            x0_time = gamma_prev_time;
            x0_exp = gamma_prev_exp;
            x0_step = gamma_prev_step;
        end
        
        % Solve nonlinear system of equations
        gamma_prev_fixed = gamma_curr_fixed;
        gamma_prev_time = gamma_curr_time;
        gamma_prev_exp = gamma_curr_exp;
        gamma_prev_step = gamma_curr_step;
        
        [gamma_curr_fixed,~,~,~] = fsolve(fun,x0_fixed,options);
        [gamma_curr_time,~,~,~] = fsolve(fun,x0_time,options);
        [gamma_curr_exp,~,~,~] = fsolve(fun,x0_exp,options);
        [gamma_curr_step,~,~,~] = fsolve(fun,x0_step,options);
               
        % Calculate actual probabilities of acceptance of passengers using their true parameters and the above price

        % Under assumption that true passenger follows assumed CPT model (same functional form) but with different parameters
        p_sR_prev_fixed = p_sR_curr_fixed;
        p_sR_prev_time = p_sR_curr_time;
        p_sR_prev_exp = p_sR_curr_exp;
        p_sR_prev_step = p_sR_curr_step;
        
        p_sR_curr_fixed = p_sR1(gamma_curr_fixed,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        p_sR_curr_time = p_sR1(gamma_curr_time,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        p_sR_curr_exp = p_sR1(gamma_curr_exp,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        p_sR_curr_step = p_sR1(gamma_curr_step,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        
        % p_sR values calculated via CPT above will differ for each trip/ride offer being considered
        % Get single value by averaging across all the samples (for convergence plot)
        p_sR_mean_fixed(i) = mean(p_sR_curr_fixed);
        p_sR_mean_time(i) = mean(p_sR_curr_time);
        p_sR_mean_exp(i) = mean(p_sR_curr_exp);
        p_sR_mean_step(i) = mean(p_sR_curr_step);
        
        error_mean_fixed(i) = abs(p_sR_mean_fixed(i) - p_star(i));
        error_mean_time(i) = abs(p_sR_mean_time(i) - p_star(i));
        error_mean_exp(i) = abs(p_sR_mean_exp(i) - p_star(i));
        error_mean_step(i) = abs(p_sR_mean_step(i) - p_star(i));
                        
        error_norm_fixed(i) = norm(p_sR_curr_fixed - p_star_vec);
        error_norm_time(i) = norm(p_sR_curr_time - p_star_vec);
        error_norm_exp(i) = norm(p_sR_curr_exp - p_star_vec);
        error_norm_step(i) = norm(p_sR_curr_step - p_star_vec);
        
        n = n+1;
        
    % After initial condition, update tariff using only feedback until p_star changes
    else
        % Gradient descent update (output feedback signal to modify input)
        
        % Learning rate schedules 
        % Time-based decay
        lr_time = lr_fixed/(1 + time_factor*n);

        % Exponential decay
        lr_exp = lr_fixed*exp(-exp_factor*n);
        
        % Step decay
        lr_step = lr_fixed*step_factor^(floor(n/step_iters));

        % Gradient of error w.r.t tariff (non-normalized)
        grad_fixed = (p_sR_curr_fixed - p_sR_prev_fixed)./(gamma_curr_fixed - gamma_prev_fixed);
        grad_time = (p_sR_curr_time - p_sR_prev_time)./(gamma_curr_time - gamma_prev_time);
        grad_exp = (p_sR_curr_exp - p_sR_prev_exp)./(gamma_curr_exp - gamma_prev_exp);
        grad_step = (p_sR_curr_step - p_sR_prev_step)./(gamma_curr_step - gamma_prev_step);

        % Set elements (trips) with gamma_currS = gamma_prevS to have zero graident (i.e. no update made)
        grad_fixed(isnan(grad_fixed)) = 0;
        grad_time(isnan(grad_time)) = 0;
        grad_exp(isnan(grad_exp)) = 0;
        grad_step(isnan(grad_step)) = 0;
                
        gamma_prev_fixed = gamma_curr_fixed;
        gamma_prev_time = gamma_curr_time;
        gamma_prev_exp = gamma_curr_exp;
        gamma_prev_step = gamma_curr_step;
        
        update_fixed = (p_sR_curr_fixed - p_star_vec).*grad_fixed;
        update_time = (p_sR_curr_time- p_star_vec).*grad_time;
        update_exp = (p_sR_curr_exp - p_star_vec).*grad_exp;
        update_step = (p_sR_curr_step - p_star_vec).*grad_step;
        
        % Norm clipping
        update_fixed = update_fixed./norm(update_fixed);
        update_time = update_time./norm(update_time);
        update_exp = update_exp./norm(update_exp);
        update_step = update_step./norm(update_step);
        
        gamma_curr_fixed = gamma_curr_fixed - lr_fixed.*update_fixed;
        gamma_curr_time = gamma_curr_time - lr_time.*update_time;
        gamma_curr_exp = gamma_curr_exp - lr_exp.*update_exp;   
        gamma_curr_step = gamma_curr_step - lr_step.*update_step;  
        
        p_sR_curr_fixed = p_sR1(gamma_curr_fixed,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        p_sR_curr_time = p_sR1(gamma_curr_time,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        p_sR_curr_exp = p_sR1(gamma_curr_exp,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        p_sR_curr_step = p_sR1(gamma_curr_step,alpha_plus_true,alpha_minus_true,beta_plus_true,beta_minus_true,lambda_true);
        
        p_sR_mean_fixed(i) = mean(p_sR_curr_fixed);
        p_sR_mean_time(i) = mean(p_sR_curr_time);
        p_sR_mean_exp(i) = mean(p_sR_curr_exp);
        p_sR_mean_step(i) = mean(p_sR_curr_step);
       
        error_mean_fixed(i) = abs(p_sR_mean_fixed(i) - p_star(i));
        error_mean_time(i) = abs(p_sR_mean_time(i) - p_star(i));
        error_mean_exp(i) = abs(p_sR_mean_exp(i) - p_star(i));
        error_mean_step(i) = abs(p_sR_mean_step(i) - p_star(i));
        
        error_norm_fixed(i) = norm(p_sR_curr_fixed - p_star_vec);
        error_norm_time(i) = norm(p_sR_curr_time - p_star_vec);
        error_norm_exp(i) = norm(p_sR_curr_exp - p_star_vec);
        error_norm_step(i) = norm(p_sR_curr_step - p_star_vec);
        
        n = n + 1;
    end
         
end
toc

% Plots
fprintf('New_gamma')
figure(1)
hold on
plot(t,p_star,'k-','LineWidth',1.2);
plot(t,p_sR_mean_fixed,'b:','LineWidth',1.2);
plot(t,p_sR_mean_time,'g-.','LineWidth',1.2);
plot(t,p_sR_mean_exp,'r--','LineWidth',1.2);
plot(t,p_sR_mean_step,'k--','LineWidth',1.2);
ylim([0 1]);
xlabel('Time (s)','Interpreter','latex');
ylabel('Probability of acceptance','Interpreter','latex');
legend('Desired $p^*(t)$','Mean $p^s_R(t)$: Fixed $\nu$','Mean $p^s_R(t)$: Time-decayed $\nu$','Mean $p^s_R(t)$: Exponentially-decayed $\nu$','Mean $p^s_R(t)$: Step-decayed $\nu$','Interpreter','latex');

figure(2)
hold on
plot(t,error_mean_fixed,'b:','LineWidth',1.2);
plot(t,error_mean_time,'g-.','LineWidth',1.2);
plot(t,error_mean_exp,'r--','LineWidth',1.2);
plot(t,error_mean_step,'k-','LineWidth',1.2);
legend('Fixed $\nu$','Time-decayed $\nu$','Exponentially-decayed $\nu$','Step-decayed $\nu$','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('Mean $|\bar{p}^s_R - \bar{p}^*|$','Interpreter','latex');

figure(3)
hold on
plot(t,error_norm_fixed,'b:','LineWidth',1.2);
plot(t,error_norm_time,'g-.','LineWidth',1.2);
plot(t,error_norm_exp,'r--','LineWidth',1.2);
plot(t,error_norm_step,'k-','LineWidth',1.2);
legend('Fixed $\nu$','Time-decayed $\nu$','Exponentially-decayed $\nu$','Step-decayed $\nu$','Interpreter','latex');
xlabel('Time (s)','Interpreter','latex');
ylabel('$\|\textbf{p}^s_R - \textbf{p}^*\|_2$','Interpreter','latex');

