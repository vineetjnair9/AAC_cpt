%% Financial risk parameters (from lottery questions)
clear;
T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_6-10_n955');
table_orig = T.Fulllaunch610n955;

% Vary if table columns change
table = convertvars(table_orig,[7:8,10:12,48:139],'double');
num_responses = size(table,1); % Total no. of survey responses

% Each row gives risk attitude parameters - alpha, beta+, beta-, lambda
cpt_lottery = zeros(num_responses,4);
error = zeros(num_responses,1); % Squared 2-norm of residual

lb = [0,0,0,1];
ub = [1,1,1,100];
x0 = [0.5,0.5,0.5,1.5];
step = [0.05,0.05,0.05,1];

options = optimoptions('lsqnonlin','Display','off');
obj_fun = @(x,respondent_num) lottery_obj(x,respondent_num,table);

%%
tic
for i = 1:50
    
    objective = @(x) obj_fun(x,i);
    
    [cpt_lottery(i,:),error(i)] = lsqnonlin(objective,x0,lb,ub,options);
   
    % Using grid search
%     [cpt_lottery(i,:),error(i)] = gs(objective,lb,ub,step);
    
    % Identifying erroneous responses

    % Testing for reflection effect

    % Testing for probability weighting
end
toc
        
%% Scatter plots
% respondents = 1:1:size(table,1);

num = 50;
respondents = 1:1:num;
cpt_lottery = cpt_lottery(1:num,:);
error = error(1:num,:);

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

%% Histograms
figure(1)
histogram(cpt_lottery(:,1),10);
xlabel('\alpha')
ylabel('Number of respondents');

figure(2)
histogram(cpt_lottery(:,2),10);
xlabel('\beta^+')
ylabel('Number of respondents');

figure(3)
histogram(cpt_lottery(:,3),10);
xlabel('\beta^-')
ylabel('Number of respondents');

figure(4)
histogram(cpt(:,4),49);
xlabel('\lambda')
ylabel('Number of respondents');

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
