%% Financial risk parameters (from lottery questions)
% Testing hypothesis 2: Approach of optimizing R for each passenger

clear; clc;

T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_6-10_n955');
table_orig = T.Fulllaunch610n955;

table = convertvars(table_orig,[7:8,10:12,48:139],'double');
num_responses = size(table,1); % Total no. of survey responses
load('valid_indices.mat');

%% Estimating CPT parameters only using valid responses
% Each row gives risk attitude parameters
num_valid = length(valid_indices);
cpt = zeros(num_valid,6);
error = zeros(num_valid,1); % Squared 2-norm of residual

lb = [0,0,0,0,0,-0.5];
ub = [1,1,1,1,1,0.5];
x0 = [0.5,0.5,0.5,0.5,0.5,0];

options = optimoptions('lsqnonlin','Display','off','FiniteDifferenceStepSize',1e-2,'DiffMinChange',1e-2,'DiffMaxChange',0.1,...
    'FiniteDifferenceType','central','MaxFunctionEvaluations',5000,'MaxIterations',4000);
weight_type = 'original';
cdf = true;
obj_fun = @(x,respondent_num) lottery_obj(x,respondent_num,table,weight_type,cdf);

tic
for i = 1:length(valid_indices)   
    j = valid_indices(i);
    objective = @(x) obj_fun(x,j);
    [cpt(i,:),error(i)] = lsqnonlin(objective,x0,lb,ub,options); 
end
toc
cpt(:,5) = cpt(:,5) * 100;
cpt(:,6) = cpt(:,6) * 10;

%% Objective function for each combination of parameter values under consideration
function F = lottery_obj(x,respondent_num,table,weight_type,cdf)
%     alpha_plus = x(1); alpha_minus = x(2); beta_plus = x(3); beta_minus = x(4); lambda = x(5); R = x(6)

    % Normalize lambda so that all 4 parameters lie be tween 0 and 1
    x(5) = x(5) * 100;
    x(6) = x(6) * 100;
    
    R = x(6);
    x = x(1:5);
    
    error = @(x,CE_actual,u1,u2,p,R) calc_error(x,CE_actual,u1,u2,p,R,weight_type,cdf);
    
    i = respondent_num;
    
    % Lottery 1
    CE = table.lottery1(i);
    u1 = 10; u2 = 100; p = 0.1;
    CE_actual = value_func(x(3),x(4),x(5),CE,R);
    f1 = error(x,CE_actual,u1,u2,p,R);
    
    % Lottery 2
    CE = table.lottery2(i);
    u1 = 0; u2 = 100; p = 0.4;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f2 = error(x,CE_actual,u1,u2,p,R);

    % Lottery 3
    CE = table.lottery3(i);
    u1 = 0; u2 = 100; p = 0.1;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f3 = error(x,CE_actual,u1,u2,p,R);

    % Lottery 4
    CE = table.lottery4(i);
    u1 = 0; u2 = 10000; p = 0.4;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f4 = error(x,CE_actual,u1,u2,p,R);

    % Lottery 5
    CE = table.lottery5(i);
    u1 = 0; u2 = 100; p = 0.9;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f5 = error(x,CE_actual,u1,u2,p,R);
    
    % Lottery 6
    CE = table.lottery6(i);
    u1 = 0; u2 = 400; p = 0.4;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f6 = error(x,CE_actual,u1,u2,p,R);

    % Lottery 7
    CE = -table.lottery7(i);
    u1 = -80; u2 = 0; p = 0.6;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f7 = error(x,CE_actual,u1,u2,p,R);

    % Lottery 8
    CE = -table.lottery8(i);
    u1 = -100; u2 = 0; p = 0.6;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f8 = error(x,CE_actual,u1,u2,p,R);

    % Lottery 9
    CE = 0;
    u1 = -25; u2 = table.lottery9(i); p = 0.5;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f9 = error(x,CE_actual,u1,u2,p,R);

    % Lottery 10
    CE = 0;
    u1 = -100; u2 = table.lottery10(i); p = 0.5;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f10 = error(x,CE_actual,u1,u2,p,R);
    
    F = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10];

end

%%
function error = calc_error(x,CE_actual,u1,u2,p,R,weight_type,cdf)
    alpha_plus = x(1); alpha_minus = x(2); beta_plus = x(3); beta_minus = x(4); lambda = x(5);
    CE_pred = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p,R,weight_type,cdf);
    
    % With normalization
    if (CE_actual == 0)
        error = abs(CE_pred - CE_actual);
    else         
        error = abs((CE_pred - CE_actual)/max(abs(u1),abs(u2)));
    end
end

%% Calculating subjective utilities
function U_sR = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p1,R,weight_type,cdf) 

    V = @(u) value_func(beta_plus,beta_minus,lambda,u,R);
    u = [u1,u2];
    p_vals = [p1,1-p1];

    [u_low,ind_low] = min(u);
    [u_high,ind_high] = max(u);

    p = p_vals(ind_low);

    w = weights(p,u_low,u_high,R,alpha_plus,alpha_minus,weight_type,cdf);
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
function w = weights(p,u1,u2,R,alpha_plus,alpha_minus,weight_type,cdf)
    F = @(u) bernoulli_cdf(p,u,u1,u2);
    pi_gain = @(x) distort_gain(x,alpha_plus,weight_type);
    pi_loss = @(x) distort_loss(x,alpha_minus,weight_type);
    
    if (cdf)
        if (u1 < R)
            w(1) = pi_loss(F(u1));
        else
            w(1) = 1 - pi_gain(1-F(u1));
        end

        if (u2 < R) 
            w(2) = pi_loss(F(u2)) - pi_loss(F(u1));
        else
            w(2) = pi_gain(1-F(u1)) - pi_gain(1-F(u2));
        end
    else
        if (u1 < R)
            w(1) = pi_loss(p);
        else
            w(1) = pi_gain(p);
        end
    
        if (u2 < R)
            w(2) = pi_loss(1-p);
        else
            w(2) = pi_gain(1-p);
        end
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
        
%% Probability weighting function - different in loss and gain regimes
function pi = distort_gain(p,alpha_plus,weight_type)
    if strcmp(weight_type,'original')
        pi = exp(-(-log(p))^alpha_plus);
    elseif strcmp(weight_type,'modified')
        pi = (p^alpha_plus)/(p^alpha_plus + (1-p)^alpha_plus)^(1/alpha_plus);
    end
end

function pi = distort_loss(p,alpha_minus,weight_type)
    if strcmp(weight_type,'original')
        pi = exp(-(-log(p))^alpha_minus);
    elseif strcmp(weight_type,'modified')
        pi = (p^alpha_minus)/(p^alpha_minus + (1-p)^alpha_minus)^(1/alpha_minus);
    end
end
            