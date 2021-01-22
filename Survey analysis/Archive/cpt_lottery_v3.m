%% Financial risk parameters (from lottery questions)
% Problem-based approach for solving nonlinear system of equations

% But optimproblem runs into an error since we're using optimization variables in exponents and if-statements

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

obj_fun = @(x,respondent_num) lottery_obj(x,respondent_num,table);

%%
x = optimvar('x',4,"LowerBound",lb,"UpperBound",ub);
for i = 1:10
    [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10] = obj_fun(x,i);
    eqn1 = f1 == 0;
    eqn2 = f2 == 0;
    eqn3 = f3 == 0;
    eqn4 = f4 == 0;
    eqn5 = f5 == 0;
    eqn6 = f6 == 0;
    eqn7 = f7 == 0;
    eqn8 = f8 == 0;
    eqn9 = f9 == 0;
    eqn10 = f10 == 0;
    
    prob = eqnproblem;
    prob.Equations.eqn1 = eqn1;
    prob.Equations.eqn2 = eqn2;
    prob.Equations.eqn3 = eqn3;
    prob.Equations.eqn4 = eqn4;
    prob.Equations.eqn5 = eqn5;
    prob.Equations.eqn6 = eqn6;
    prob.Equations.eqn7 = eqn7;
    prob.Equations.eqn8 = eqn8;
    prob.Equations.eqn9 = eqn9;
    prob.Equations.eqn10 = eqn10;
    x0.x = [0.5,0.5,0.5,1.5];
    
    cpt_lottery(i,:) = solve(prob,x0);
end
        
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

%% Objective function for each combination of parameter values under consideration
function [f1,f2,f3,f4,f5,f6,f7,f8,f9,f10] = lottery_obj(x,respondent_num,table)
%     alpha = x(1); beta_plus = x(2); beta_minus = x(3); lambda = x(4);

    R = 0;
    i = respondent_num;
    % Lottery 1
    CE = table.lottery1(i);
    u1 = 10; u2 = 100; p = 0.1;
    CE_sR = value_func(x(2),x(3),x(4),CE,R);
    f1 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);
    
    % Lottery 2
    CE = table.lottery2(i);
    u1 = 0; u2 = 100; p = 0.4;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f2 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);

    % Lottery 3
    CE = table.lottery3(i);
    u1 = 0; u2 = 100; p = 0.1;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f3 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);

    % Lottery 4
    CE = table.lottery4(i);
    u1 = 0; u2 = 10000; p = 0.4;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f4 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);

    % Lottery 5
    CE = table.lottery5(i);
    u1 = 0; u2 = 100; p = 0.9;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f5 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);
    
    % Lottery 6
    CE = table.lottery6(i);
    u1 = 0; u2 = 400; p = 0.4;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f6 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);

    % Lottery 7
    CE = table.lottery7(i);
    u1 = -80; u2 = 0; p = 0.6;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f7 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);

    % Lottery 8
    CE = table.lottery8(i);
    u1 = -100; u2 = 0; p = 0.6;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f8 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);

    % Lottery 9
    CE = 0;
    u1 = -25; u2 = table.lottery9(i); p = 0.5;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f9 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);

    % Lottery 10
    CE = 0;
    u1 = -100; u2 = table.lottery10(i); p = 0.5;
    CE_sR = -value_func(x(2),x(3),x(4),CE,R);
    f10 = calc_error(x(1),x(2),x(3),x(4),CE_sR,u1,u2,p,R);
end

%%
function error = calc_error(alpha,beta_plus,beta_minus,lambda,CE_sR,u1,u2,p,R)
    CE_pred = subjective(alpha,beta_plus,beta_minus,lambda,u1,u2,p,R,'original');
    error = CE_pred - CE_sR;
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
        V = (u-R).^beta_plus;
    else
        V = -lambda*(R-u).^beta_minus;
    end
end

%% Calculating subjective weights
function w = weights(p,u1,u2,R,alpha)
    F = @(u) bernoulli_cdf(p,u,u1,u2);
    pi = @(p) distort(p,alpha);
    
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
function pi = distort(p,alpha)
    pi = exp(-(-log(p))^alpha);
end
