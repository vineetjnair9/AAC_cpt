%% CPT parameters derivation

% Changes from 1st version
% (1) Using modified approach of checking verifying whether it's a gain or loss
% and solving 4x4 non linear system for all 4 parameters at once
% (2) Using self-reference with SMoDS itself rather than current most common mode
% (3) Simplified calculation of probability weights

% paramhat = load('paramhat.mat');
% paramhat = paramhat.paramhat;

%% Mode-choice parameters

walk_time = paramhat(1);
wait_time = paramhat(2);
drive_time = paramhat(3);
transit_time = paramhat(4);
rideshare_time = paramhat(5);
price = paramhat(6);
ASC_drive = paramhat(7);
ASC_transit = paramhat(8);
ASC_exclusive = paramhat(9);
ASC_pooled = paramhat(10);

% Value of time spent on different modes (as % of hourly wage)
% Need to adjust regression to account for wages (normalize price variable by income/hourly wage)
% Annual income --> Hourly wage

VOT_walking = abs(price/walk_time)*100
VOT_waiting = abs(price/wait_time)*100
VOT_drive = abs(price/drive_time)*100
VOT_transit = abs(price/transit_time)*100
VOT_rideshare = abs(price/rideshare_time)*100

%% Risk scenarios

T = load('Pilotdata0117');
table_orig = T.Pilotdata117;

% Vary if table columns change
table = convertvars(table_orig,[2:5,19:26,32:116],'double');

%%

num_responses = size(table,1); % Total no. of survey responses
scenarios = 5*num_responses; % 5 choice scenarios per respondent
alternatives = scenarios*4; % 4 alternatives per choice situation
X = zeros(alternatives,13); % Observations for each scenario on 4 predictor variables

% Each row gives risk attitude parameters - alpha, beta+, beta-, lambda
cpt = zeros(num_responses,4);

for i = 1:num_responses
    
    % Calculate static self-reference for each passenger based on average price & travel time 
    SMODS_cost_base = 2.2 + table.Pool_cost(i)*table.Distance(i);
    R = walk_time * 3.5 + wait_time * 3.5 + drive_time * 8 * table.Distance(i) + price * SMODS_cost_base + ASC_pooled;
    
    % Reference = Driving
    if table.Reference(i) == 'Driving'
        drive_cost = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Distance(i);
              
        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * 5 + wait_time * 0 + drive_time * table.Driving_ref_1(i) + price * drive_cost + ASC_drive;
        
        SMODS_cost = (2.2 + table.Pool_cost(i)*table.Distance(i))*0.8; 
        u1 = walk_time * 2 + wait_time * 0 + rideshare_time * 7 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.1;

        u2 = walk_time * 3 + wait_time * 1 + rideshare_time * 6 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.9;  
        
        f1 = @(alpha,beta,gamma,lambda) subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R) - CE;
        
        % Scenario 2
        CE = walk_time * 5 + wait_time * 0 + drive_time * table.Driving_ref_2(i) + price * drive_cost + ASC_drive;
        
        p1 = 0.6;
        p2 = 0.4;
        
        f2 = @(alpha,beta,gamma,lambda) subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R) - CE;
                      
        % Scenario 3
        CE = walk_time * 0 + wait_time * 0 + drive_time * table.Driving_ref_3(i) + price * drive_cost + ASC_drive;
        
        SMODS_cost = (2.2 + table.Pool_cost(i)*table.Distance(i))*1.2; 
        
        u1 = walk_time * 7 + wait_time * 4 + rideshare_time * 6 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.8;

        u2 = walk_time * 5 + wait_time * 6 + rideshare_time * 10 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.2;
       
        f3 = @(alpha,beta,gamma,lambda) subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R) - CE;
        
        % Scenario 4
        CE = walk_time * 5 + wait_time * 0 + drive_time * table.Driving_ref_4(i) + price * drive_cost + ASC_drive;
        
        SMODS_cost = 2.2 + table.Pool_cost(i)*table.Distance(i); 

        u1 = walk_time * 7 + wait_time * 6 + rideshare_time * 9 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.5; % loss

        u2 = walk_time * 1 + wait_time * 1 + rideshare_time * 7 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.5; % gain
        
        f4 = @(alpha,beta,gamma,lambda) subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R) - CE;
        
        fun = @(x) system(f1,f2,f3,f4,x);
        
        lb = [0.01,0.01,0.01,1.01];
        ub = [0.99,0.99,0.99,inf];
        options = optimoptions(@lsqnonlin,'OptimalityTolerance',1e-10,'FiniteDifferenceStepSize',1e-2);
        cpt(i,:) = lsqnonlin(fun,[0.5,0.5,0.5,1.5],lb,ub,options);

%         cpt(i,:) = fsolve(fun,[0.5,0.5,0.5,1.5])
    elseif (table.Reference(i) == 'Public transit' && table.transit_mode(i) == 'Subway (T)')
        Subway_cost = 2.65;
        
        transit_walk = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * 1.05;
        
        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * transit_walk + wait_time * 15 + transit_time * table.Subway_ref_1(i) + price * Subway_cost + ASC_transit;
        
        SMODS_cost = (2.2 + table.Pool_cost(i)*table.Distance(i))*0.8; 
        u1 = walk_time * 2 + wait_time * 0 + rideshare_time * 6 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.1;

        u2 = walk_time * 3 + wait_time * 1 + rideshare_time * 7 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.9;  
        
        f1 = @(alpha,beta,gamma,lambda) subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R) - CE;
        
        % Scenario 2
        CE = walk_time * transit_walk + wait_time * 15 + transit_time * table.Subway_ref_2(i) + price * Subway_cost + ASC_transit;
        
        p1 = 0.6;
        p2 = 0.4;
        
        f2 = @(alpha,beta,gamma,lambda) subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R) - CE;
                      
        % Scenario 3
        transit_walk = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * 0.95;
        
        CE = walk_time * transit_walk + wait_time * 2 + transit_time * table.Subway_ref_3(i) + price * Subway_cost + ASC_transit;
        
        SMODS_cost = (2.2 + table.Pool_cost(i)*table.Distance(i))*1.2; 
        
        u1 = walk_time * 7 + wait_time * 6 + rideshare_time * 6 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.8;

        u2 = walk_time * 5 + wait_time * 4 + rideshare_time * 9 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.2;
       
        f3 = @(alpha,beta,gamma,lambda) subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R) - CE;
        
        % Scenario 4
        transit_walk = (table.Walk_hometransit(i) + table.Walk_transitwork(i));
        
        CE = walk_time * transit_walk + wait_time * 2 + transit_time * table.Subway_ref_4(i) + price * Subway_cost + ASC_transit;
        
        SMODS_cost = 2.2 + table.Pool_cost(i)*table.Distance(i); 

        u1 = walk_time * 7 + wait_time * 6 + rideshare_time * 8 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.5; % loss

        u2 = walk_time * 2 + wait_time * 2 + rideshare_time * 6 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.5; % gain
        
        f4 = @(alpha,beta,gamma,lambda) subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R) - CE;
        
        fun = @(x) system(f1,f2,f3,f4,x);
        
        lb = [0.01,0.01,0.01,1.01];
        ub = [0.99,0.99,0.99,inf];
        options = optimoptions(@lsqnonlin,'OptimalityTolerance',1e-10,'FiniteDifferenceStepSize',1e-2);
        cpt(i,:) = lsqnonlin(fun,[0.5,0.5,0.5,1.5],lb,ub,options);
        
        % Also need to add in reference cases for transit-bus, transit-commuter rail & exclusive rideshare
 
    end
end        
        
%% Plotting function

%% 
function F = system(f1,f2,f3,f4,x)
% Here x(1) = alpha, x(2) = beta, x(3) = gamma, x(4) = lambda
F(1) = f1(x(1),x(2),x(3),x(4));
F(2) = f2(x(1),x(2),x(3),x(4));
F(3) = f3(x(1),x(2),x(3),x(4));
F(4) = f4(x(1),x(2),x(3),x(4));
end

%% Using new defn of weights
function U = subjective(alpha,beta,gamma,lambda,u1,u2,p1,p2,R)
u = [u1,u2];
p_vals = [p1,p2];
[u_low,ind_low] = min(u);
[u_high,ind_high] = max(u);
p = p_vals(ind_high); % prob of higher outcome
if (u_high < R) % Pure loss
    w = exp(-(-log(p))^alpha);
    U = -w .* lambda .* (R - u_low).^gamma - (1-w) .* lambda .* (R - u_high).^gamma;
elseif (u_low >= R) % Pure gain
    w = -exp(-(-log(p))^alpha);
    U = (1-w) .* (u_low - R).^beta + w .* (u_high - R).^beta;
else % Mixed prospects: u_low < R & u_high > R
    w1 = exp(-(-log(p))^alpha); % gain
    w2 = exp(-(-log(1-p))^alpha); % loss
    U = -w2 .* lambda .* (R - u_low).^gamma + w1 .* (u_high - R).^beta;
end
end

%% Using original defn of weights
function U = sub_util(alpha,beta,gamma,lambda,u1,u2,R,pd)
u_low = min(u1,u2);
u_high = max(u1,u2);
if (u_high < R) % Pure loss
    w1 = exp(-(-log(cdf(pd,u_low)))^alpha);
    w2 = exp(-(-log(cdf(pd,u_high)))^alpha) - exp(-(-log(cdf(pd,u_low)))^alpha);
    U = -w1 .* lambda .* abs(R - u_low).^gamma - w2 .* lambda .* abs(R - u_high).^gamma - R;
elseif (u_low >= R) % Pure gain
    w1 = -exp(-(-log(1-cdf(pd,u_low)))^alpha);
    w2 = exp(-(-log(1-cdf(pd,u_low)))^alpha) - exp(-(-log(1-cdf(pd,u_high)))^alpha);
    U = w1 .* abs(u_low - R).^beta + w2 .* abs(u_high - R).^beta - R;
else % Mixed prospects
    w1 = exp(-(-log(cdf(pd,u_low)))^alpha); % loss
    w2 = exp(-(-log(1-cdf(pd,u_low)))^alpha) - exp(-(-log(1-cdf(pd,u_high)))^alpha); % gain
    F = -w1 .* lambda .* abs(R - u_low).^gamma + w2 .* abs(u_high - R).^beta - R;
end
end