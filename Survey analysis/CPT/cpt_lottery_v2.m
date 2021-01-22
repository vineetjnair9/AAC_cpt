%% Financial risk parameters (from lottery questions)
clear;
T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_6-10_n955');
table_orig = T.Fulllaunch610n955;
load('valid_indices.mat');

% Vary if table columns change
table = convertvars(table_orig,[7:8,10:12,48:139],'double');
num_responses = size(table,1); % Total no. of survey responses
num_valid = length(valid_indices);

% Each row gives risk attitude parameters - alpha, beta+, beta-, lambda
cpt_lottery = zeros(num_valid,4);
error = zeros(num_valid,1); % Squared 2-norm of residual

lb = [0,0,0,1];
ub = [1,1,1,100];
x0 = [0.5,0.5,0.5,1.5];
step = [0.05,0.05,0.05,1];

options = optimoptions('lsqnonlin','Display','off');
obj_fun = @(x,respondent_num) lottery_obj(x,respondent_num,table);

%%
tic
parfor i = 1:length(valid_indices)
    j = valid_indices(i);
    objective = @(x) obj_fun(x,j);
    
    [cpt_lottery(i,:),error(i)] = lsqnonlin(objective,x0,lb,ub,options);
   
    % Using grid search
%     [cpt_lottery(i,:),error(i)] = gs(objective,lb,ub,step);

end
toc
        
%% Scatter plots

num = num_valid;
respondents = 1:1:num;

figure(1)
scatter(respondents,cpt_lottery(:,1))
ylabel('$\alpha$','Interpreter','latex')
xlabel('Respondent','Interpreter','latex')

figure(2)
scatter(respondents,cpt_lottery(:,2))
ylabel('$\beta^+$','Interpreter','latex')
xlabel('Respondent','Interpreter','latex')

figure(3)
scatter(respondents,cpt_lottery(:,3))
ylabel('$\beta^-$','Interpreter','latex')
xlabel('Respondent','Interpreter','latex')

figure(4)
scatter(respondents,cpt_lottery(:,4))
ylabel('$\lambda$','Interpreter','latex')
xlabel('Respondent','Interpreter','latex')

figure(5)
scatter(respondents,error)
ylabel('Sum of squared errors','Interpreter','latex')
xlabel('Respondent','Interpreter','latex')

%% Histograms
figure(1)
histogram(cpt_lottery(:,1));
xlabel('$\alpha$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')

figure(2)
histogram(cpt_lottery(:,2));
xlabel('$\beta^+$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')

figure(3)
histogram(cpt_lottery(:,3));
xlabel('$\beta^-$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')

figure(4)
histogram(cpt_lottery(:,4));
xlabel('$\lambda$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')

%% Objective function for each combination of parameter values under consideration
function F = lottery_obj(x,respondent_num,table)
%     alpha = x(1); beta_plus = x(2); beta_minus = x(3); lambda = x(4);

    R = 0;
    i = respondent_num;
    
    % Lottery 1
    CE = table.lottery1(i);
    u1 = 10; u2 = 100; p = 0.1;
%     CE_actual = CE;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f1 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);
    
    % Lottery 2
    CE = table.lottery2(i);
    u1 = 0; u2 = 100; p = 0.4;
%     CE_actual = CE;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f2 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);

    % Lottery 3
    CE = table.lottery3(i);
    u1 = 0; u2 = 100; p = 0.1;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f3 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);

    % Lottery 4
    CE = table.lottery4(i);
    u1 = 0; u2 = 10000; p = 0.4;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f4 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);

    % Lottery 5
    CE = table.lottery5(i);
    u1 = 0; u2 = 100; p = 0.9;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f5 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);
    
    % Lottery 6
    CE = table.lottery6(i);
    u1 = 0; u2 = 400; p = 0.4;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f6 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);

    % Lottery 7
    CE = -table.lottery7(i);
    u1 = -80; u2 = 0; p = 0.6;
%     CE_actual = CE;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f7 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);

    % Lottery 8
    CE = -table.lottery8(i);
    u1 = -100; u2 = 0; p = 0.6;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f8 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);

    % Lottery 9
    CE = 0;
    u1 = -25; u2 = table.lottery9(i); p = 0.5;
%     CE_actual = CE;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f9 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);

    % Lottery 10
    CE = 0;
    u1 = -100; u2 = table.lottery10(i); p = 0.5;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f10 = calc_error(x(1),x(2),x(3),x(4),CE_actual,u1,u2,p,R);
    
    % Using all 10 lotteries
    F = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10];

    % Using only 4 lottery equations (in 4 unknowns)
%     F = [f1; f2; f7; f9];

end

%%
function error = calc_error(alpha,beta_plus,beta_minus,lambda,CE_actual,u1,u2,p,R)
    CE_pred = subjective(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R,'original');
    
    % No normalization 
    error = CE_pred - CE_actual;
    
%     % With normalization
%     if (CE_actual == 0)
%         error = abs(CE_pred - CE_actual);
%     else         
%         % Taking relative error w.r.t true value (CE_sR)
%         error = abs((CE_pred - CE_sR)/CE_sR);
%     end
end

%% Calculating subjective utilities
function U_sR = subjective(alpha,beta_plus,beta_minus,lambda,u1,u2,p1,R,weight_type) 

    V = @(u) value_func(beta_plus,beta_minus,lambda,u,R);
    u = [u1,u2];
    p_vals = [p1,1-p1];

    [u_low,ind_low] = min(u);
    [u_high,ind_high] = max(u);

    p = p_vals(ind_low);

    w = weights(p,u_low,u_high,R,alpha,weight_type);
    U_sR = w(1) * V(u_low) + w(2) * V(u_high);
end
    
%% Value function
function V = value_func(beta_plus,beta_minus,lambda,u,R)
    if (u >= R)
        V = (u-R)^beta_plus;
    else
        V = -lambda*(R-u)^beta_minus;
    end
end

%% Calculating subjective weights
function w = weights(p,u1,u2,R,alpha,weight_type)
    F = @(u) bernoulli_cdf(p,u,u1,u2);
    pi = @(x) distort(x,alpha,weight_type);
    
    if (u1 < R)
        w(1) = pi(F(u1));
    else
        w(1) = 1 - pi(1-F(u1));
    end
    
    if (u2 < R) 
        w(2) = pi(F(u2)) - pi(F(u1));
    else
        w(2) = pi(1-F(u1)) - pi(1-F(u2));
    end
end

%% Bernoulli CDF (binomial distribution with only 1 trial & 2 outcomes)
function F = bernoulli_cdf(p,u,u1,u2)
    if (u < u1)
        F = 0;
    elseif (u >= u1 && u < u2)
        F = p;
    else
        F = 1;
    end
end
        
%% Probability weighting function
function pi = distort(p,alpha,weight_type)
    if strcmp(weight_type,'original')
        pi = exp(-(-log(p))^alpha);
    elseif strcmp(weight_type,'modified')
        pi = (p^alpha)/(p^alpha + (1-p)^alpha)^(1/alpha);
    end
end

%% Grid search function
function [opt_params,resnorm] = gs(fun,lb,ub,step)
    alpha_vals = lb(1):step(1):ub(1);
    beta_plus_vals = lb(2):step(2):ub(2);
    beta_minus_vals = lb(3):step(3):ub(3);
    lambda_vals = lb(4):step(4):ub(4);
    
    residuals = zeros(length(alpha_vals),length(beta_plus_vals),length(beta_minus_vals),length(lambda_vals));
    for h = 1:length(alpha_vals)
        for i = 1:length(beta_plus_vals)
            for j = 1:length(beta_minus_vals)
                for k = 1:length(lambda_vals)
                    x(1) = alpha_vals(h); x(2) = beta_plus_vals(i);
                    x(3) = beta_minus_vals(j); x(4) = lambda_vals(k);
                    residuals(h,i,j,k) = norm(fun(x));
                end
            end
        end
    end
    [resnorm, minidx] = min(residuals(:));
    [h, i, j, k] = ind2sub(size(residuals), minidx);
    opt_params = [alpha_vals(h),beta_plus_vals(i),beta_minus_vals(j),lambda_vals(k)];
end
            