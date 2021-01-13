%% Financial risk parameters (from lottery questions)
clear; clc;
cd 'C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC'
T = load('Data\Panel\Full_launch_6-10_n955');
table_orig = T.Fulllaunch610n955;

% Vary if table columns change
table = convertvars(table_orig,[7:8,10:12,48:139],'double');
num_responses = size(table,1); % Total no. of survey responses

% Each row gives risk attitude parameters
cpt_lottery = zeros(num_responses,5);
error = zeros(num_responses,1); % Squared 2-norm of residual

lb = [0,0,0,0,0];
ub = [1,1,1,1,1];
x0 = [0.5,0.5,0.5,0.5,0.5];

options = optimoptions('lsqnonlin','Display','off','FiniteDifferenceStepSize',1e-2,'DiffMinChange',1e-2,'DiffMaxChange',0.1,...
    'FiniteDifferenceType','central','MaxFunctionEvaluations',5000,'MaxIterations',4000);
obj_fun = @(x,respondent_num) lottery_obj(x,respondent_num,table);

%%
valid = zeros(num_responses,1);

for i = 1:20/num_responses
    % Identifying erroneous responses
    CE_low = table.lottery5(i); CE_med = table.lottery2(i); CE_high = table.lottery3(i); 
    
    if (CE_low <= CE_med) && (CE_med <= CE_high) && (table.lottery7(i) <= table.lottery8(i))
        valid(i) = 1;
    end
end

% Only use valid/rational responses going forward
valid_indices = find(valid);
num_valid = length(valid_indices);
reflection_effect = zeros(num_valid,1);
prob_weighting1 = zeros(num_valid,1);
prob_weighting2 = zeros(num_valid,1);
prob_weighting3 = zeros(num_valid,1);
prob_weighting_any = zeros(num_valid,1);
cpt_loss_aversion = zeros(num_valid,1);

for i = 1:length(valid_indices)   
    j = valid_indices(i);
    
    % All these tests are borrowed from Rieger et al. 2017
    % Testing for reflection effect
    if (table.lottery2(j) <= 60) && (table.lottery8(j) <= 60)
        reflection_effect(i) = 1;
    end
   
    % Testing for probability weighting
    if (table.lottery5(j)/table.lottery2(j)) > (1/6)
        prob_weighting1(i) = 1;
    end
    
    if (table.lottery2(j)/table.lottery3(j)) > (2/3)
        prob_weighting2(i) = 1;
    end
    
    if (table.lottery5(j)/table.lottery3(j)) > (1/9)
        prob_weighting3(i) = 1;
    end
    
    if (prob_weighting1(i) == 1) || (prob_weighting2(i) == 1) || (prob_weighting3(i) == 1)
        prob_weighting_any(i) = 1;
    end
    
    % Testing whether loss aversion in PT is more valid than PT
    cpt_loss_aversion(i) = mean([table.lottery9(j)/25, table.lottery10(j)/100]);
end

% Percentage of valid responses
valid_responses = (sum(valid)/num_responses)*100

% Percentage of valid respondents displaying reflection effect
reflection_responses = (sum(reflection_effect)/num_valid)*100

% Percentage of valid respondents displaying probability weighting
prob_weighting1_ratio = (sum(prob_weighting1)/num_valid)*100
prob_weighting2_ratio = (sum(prob_weighting2)/num_valid)*100
prob_weighting3_ratio = (sum(prob_weighting3)/num_valid)*100
prob_weighting_any_ratio = (sum(prob_weighting_any)/num_valid)*100

% Loss aversion (A/B) ratio
mean_ratio = mean(cpt_loss_aversion)

%% Estimating CPT parameters on all responses
tic
for i = 1:num_responses
    objective = @(x) obj_fun(x,i);
    [cpt_lottery(i,:),error(i)] = lsqnonlin(objective,x0,lb,ub,options); 
end
toc

%% Estimating CPT parameters only using valid responses
tic
for i = 1:length(valid_indices)   
    j = valid_indices(i);
    objective = @(x) obj_fun(x,j);
    [cpt_lottery(i,:),error(i)] = lsqnonlin(objective,x0,lb,ub,options); 
end
toc

%% With normalized lambda
cpt_lottery(:,5) = cpt_lottery(:,5) * 100;

%% Objective function for each combination of parameter values under consideration
function F = lottery_obj(x,respondent_num,table)
%     alpha_plus = x(1); alpha_minus = x(2); beta_plus = x(3); beta_minus = x(4); lambda = x(5);

    % Normalize lambda so that all 4 parameters lie be tween 0 and 1
    x(5) = x(5) * 100;
    
    R = 0;
    i = respondent_num;
    
    % Lottery 1
    CE = table.lottery1(i);
    u1 = 10; u2 = 100; p = 0.1;
%     CE_actual = CE;
    CE_actual = value_func(x(3),x(4),x(5),CE,R);
    f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);
    
    % Lottery 2
    CE = table.lottery2(i);
    u1 = 0; u2 = 100; p = 0.4;
%     CE_actual = CE;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);

    % Lottery 3
    CE = table.lottery3(i);
    u1 = 0; u2 = 100; p = 0.1;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);

    % Lottery 4
    CE = table.lottery4(i);
    u1 = 0; u2 = 10000; p = 0.4;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);

    % Lottery 5
    CE = table.lottery5(i);
    u1 = 0; u2 = 100; p = 0.9;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);
    
    % Lottery 6
    CE = table.lottery6(i);
    u1 = 0; u2 = 400; p = 0.4;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);

    % Lottery 7
    CE = -table.lottery7(i);
    u1 = -80; u2 = 0; p = 0.6;
%     CE_actual = CE;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f7 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);

    % Lottery 8
    CE = -table.lottery8(i);
    u1 = -100; u2 = 0; p = 0.6;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f8 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);

    % Lottery 9
    CE = 0;
    u1 = -25; u2 = table.lottery9(i); p = 0.5;
%     CE_actual = CE;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f9 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);

    % Lottery 10
    CE = 0;
    u1 = -100; u2 = table.lottery10(i); p = 0.5;
    CE_actual = value_func(x(2),x(3),x(4),CE,R);
    f10 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_actual,u1,u2,p,R);
    
    % Using all 10 lotteries
    F = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10];

    % Using only 4 lottery equations (in 4 unknowns)
%     F = [f1; f2; f7; f9];

end

%%
function error = calc_error(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,CE_actual,u1,u2,p,R)
    CE_pred = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p,R,'modified');
    
    % No normalization 
%     error = CE_pred - CE_actual;
    
    % With normalization
    if (CE_actual == 0)
        error = abs(CE_pred - CE_actual);
    else         
        % Taking relative error w.r.t true value
%         error = abs((CE_pred - CE_actual)/CE_actual);

        % Taking relative error w.r.t max outcome
        error = abs((CE_pred - CE_actual)/max(abs(u1),abs(u2)));
    end
end

%% Calculating subjective utilities
function U_sR = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p1,R,weight_type) 

    V = @(u) value_func(beta_plus,beta_minus,lambda,u,R);
    u = [u1,u2];
    p_vals = [p1,1-p1];

    [u_low,ind_low] = min(u);
    [u_high,ind_high] = max(u);

    p = p_vals(ind_low);

    w = weights(p,u_low,u_high,R,alpha_plus,alpha_minus,weight_type);
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
function w = weights(p,u1,u2,R,alpha_plus,alpha_minus,weight_type)
    F = @(u) bernoulli_cdf(p,u,u1,u2);
    pi_gain = @(x) distort_gain(x,alpha_plus,weight_type);
    pi_loss = @(x) distort_loss(x,alpha_minus,weight_type);
    
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
            