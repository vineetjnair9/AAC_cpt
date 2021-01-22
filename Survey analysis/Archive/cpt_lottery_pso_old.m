%% Financial risk parameters (from lottery questions)
clear;
T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_5-17_n307');
table_orig = T.Fulllaunch517n307;

% Vary if table columns change
table = convertvars(table_orig,[7:8,10:12,48:139],'double');

tic;
[cpt_lottery,error] = estimate_params_lottery(table,'original');
toc

%% Scatter plots
respondents = 1:1:size(table,1);

figure(1)
scatter(respondents,cpt_lottery(:,1))
ylabel('\alpha')
xlabel('Respondent')

figure(2)
scatter(respondents,cpt_lottery(:,2))
ylabel('\beta^+')
xlabel('Respondent')

figure(3)
scatter(respondents,cpt_lottery(:,3))
ylabel('\beta^-')
xlabel('Respondent')

figure(4)
scatter(respondents,cpt_lottery(:,4))
ylabel('\lambda')
xlabel('Respondent')

figure(5)
scatter(respondents,error)
ylabel('Sum of normalized errors')
xlabel('Respondent')

%%
function [cpt_lottery,error] = estimate_params_lottery(table,weight_type)
    
    num_responses = size(table,1); % Total no. of survey responses
    
    % Each row gives risk attitude parameters - alpha, beta+, beta-, lambda
    cpt_lottery = zeros(num_responses,4);
    error = zeros(num_responses,1); % Squared 2-norm of residual
   
    lb = [0,0,0,1];
    ub = [1,1,1,100];

    rng default % for reproducibility
    options = optimoptions('particleswarm','UseParallel',true);
%         options.OptimalityTolerance = 1e-14;
%         options.MaxFunctionEvaluations = Inf;
%         options.MaxIterations = Inf;
    
    R = 0;
    alpha = sym('alpha'); beta_plus = sym('beta_plus'); beta_minus = sym('beta_minus'); lambda = sym('lambda');
    
    for i = 1:5
        % Lottery 1
        CE = table.lottery1(i);
        u1 = 10; u2 = 100; p = 0.1;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f1 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f1 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));

        % Lottery 2
        CE = table.lottery2(i);
        u1 = 0; u2 = 100; p = 0.4;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f2 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f2 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));

        % Lottery 3
        CE = table.lottery3(i);
        u1 = 0; u2 = 100; p = 0.1;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f3 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f3 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));

        % Lottery 4
        CE = table.lottery4(i);
        u1 = 0; u2 = 10000; p = 0.4;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f4 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f4 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));
        
        % Lottery 5
        CE = table.lottery5(i);
        u1 = 0; u2 = 100; p = 0.9;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f5 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f5 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));
        
        % Lottery 6
        CE = table.lottery6(i);
        u1 = 0; u2 = 400; p = 0.4;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f6 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f6 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));

        % Lottery 7
        CE = table.lottery7(i);
        u1 = -80; u2 = 0; p = 0.6;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f7 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f7 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));

        % Lottery 8
        CE = table.lottery8(i);
        u1 = -100; u2 = 0; p = 0.6;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f8 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f8 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));

        % Lottery 9
        CE = 0;
        u1 = -25; u2 = table.lottery9(i); p = 0.5;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f9 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f9 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));

        % Lottery 10
        CE = 0;
        u1 = -100; u2 = table.lottery10(i); p = 0.5;
        CE_sR = @(alpha,beta_plus,beta_minus,lambda) -value_func(beta_plus,beta_minus,lambda,CE,R);
        f10 = @(alpha,beta_plus,beta_minus,lambda) calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R);
        %f10 = @(alpha,beta_plus,beta_minus,lambda) abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R) + CE_sR(alpha,beta_plus,beta_minus,lambda))/max(abs(sub(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R)), abs(CE_sR(alpha,beta_plus,beta_minus,lambda)));

        fun = @(x) system(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,x);
        [cpt_lottery(i,:),error(i)] = particleswarm(@(x) max(fun(x)),4,lb,ub,options);
        
        % Identifying erroneous responses
        
        % Testing for reflection effect
        
        % Testing for probability weighting
        
    end
end
%%
function error_normalized = calc_error(a,b,c,d,CE_sR,u1,u2,p,R)
    sub = @(a,b,c,d,u1,u2,p,R) subjective(a,b,c,d,u1,u2,p,R,'original');
    CE_pred = @(a,b,c,d) sub(a,b,c,d,u1,u2,p,R);
    if (CE_sR(a,b,c,d) == 0)
        error_normalized = abs(CE_pred(a,b,c,d) + CE_sR(a,b,c,d));
    else 
        % No normalization 
        % error_normalized = abs(CE_pred(a,b,c,d) + CE_sR(a,b,c,d));
        
        % Taking relative error with respect to largest value
        % error_normalized = abs(CE_pred(a,b,c,d) + CE_sR(a,b,c,d))/max(abs(CE_pred(a,b,c,d)),abs(CE_sR(a,b,c,d)));
    
        % Taking relative error w.r.t true value (CE_sR)
        error_normalized = abs((CE_pred(a,b,c,d) + CE_sR(a,b,c,d))/CE_sR(a,b,c,d));
    end
end

%% Set up system of nonlinear equations
function F = system(f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,x)
% Here x(1) = alpha, x(2) = beta_plus, x(3) = beta_minus, x(4) = lambda
F(1) = f1(x(1),x(2),x(3),x(4));
F(2) = f2(x(1),x(2),x(3),x(4));
F(3) = f3(x(1),x(2),x(3),x(4));
F(4) = f4(x(1),x(2),x(3),x(4));
F(5) = f5(x(1),x(2),x(3),x(4));
F(6) = f6(x(1),x(2),x(3),x(4));
F(7) = f7(x(1),x(2),x(3),x(4));
F(8) = f8(x(1),x(2),x(3),x(4));
F(9) = f9(x(1),x(2),x(3),x(4));
F(10) = f10(x(1),x(2),x(3),x(4));
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
    pi = @(u) distort(p,alpha,weight_type);
    
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
                    residuals(h,i,j,k) = abs(sum(fun(x)));
                end
            end
        end
    end
    [resnorm, minidx] = min(residuals(:));
    [h, i, j, k] = ind2sub(size(residuals), minidx);
    opt_params = [alpha_vals(h),beta_plus_vals(i),beta_minus_vals(j),lambda_vals(k)];
end
            