%% CPT parameters derivation

% Changes from 1st version
% (1) Using modified approach of checking verifying whether it's a gain or loss
% and solving 4x4 non linear system for all 4 parameters at once
% (2) Using static self-reference with SMoDS itself rather than current most common mode
% (3) Simplified calculation of probability weights

% paramhat = load('paramhat.mat');
% paramhat = paramhat.paramhat;

%% Mode-choice parameters

walk_time = paramhat(1);
wait_time = paramhat(2);
drive_time = paramhat(3);
transit_time = paramhat(4);
exc_time = paramhat(5);
pool_time = paramhat(6);
price = paramhat(7)*0.1;
ASC_drive = paramhat(8); % Or set = 0 if used as reference/baseline alternative
ASC_transit = paramhat(9);
ASC_exclusive = paramhat(10);
ASC_pooled = paramhat(11);

% Value of time (in $/min)
VOT_walking = abs(price/walk_time)
VOT_waiting = abs(price/wait_time)
VOT_drive = abs(price/drive_time)
VOT_transit = abs(price/transit_time)
VOT_exc = abs(price/exc_time)
VOT_pool = abs(price/pool_time)

%% Value of time spent on different modes (as % of hourly wage)
% Need to adjust regression to account for wages (normalize price variable by income/hourly wage)
% Annual income --> Hourly wage

VOT_walking = abs(price/walk_time)*100
VOT_waiting = abs(price/wait_time)*100
VOT_drive = abs(price/drive_time)*100
VOT_transit = abs(price/transit_time)*100
VOT_rideshare = abs(price/rideshare_time)*100

%% Risk scenarios

T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Pilot_v2_2-20-20.mat');
table_orig = T.Pilotv222020;

% Vary if table columns change
table = convertvars(table_orig,[1,4:5,7:8,43:125],'double');

%%
num_responses = size(table,1); % Total no. of survey responses

% Each row gives risk attitude parameters - alpha, beta+, beta-, lambda
cpt = zeros(num_responses,4);

for i = 1:num_responses
       
    % Calculate static self-reference for each passenger based on average price & travel time 
    SMODS_cost_base = 2.2 + table.Pool_cost(i)*table.Distance(i);
    R = walk_time * 4 + wait_time * 3.5 + pool_time * 6 * table.Distance(i) + price * SMODS_cost_base + ASC_pooled;
    
    transit_walk = table.transit_walk_1(i) + table.transit_walk_2(i);
    transit_cost = table.Transit_cost(i);
    
    p1 = table.p1(i)/100;
    p2 = table.p2(i)/100;
    p3 = table.p3(i)/100;
    p4 = table.p4(i)/100;
    
    SMODS_cost = (2.2 + table.Pool_cost(i)*table.Distance(i));
%     lb = [-inf,-inf,-inf,-inf];
%     ub = [inf,inf,inf,inf];
    
    lb = [0.1,0.1,0.1,0];
    ub = [0.99,0.99,0.99,inf];
    options = optimoptions(@lsqnonlin,'OptimalityTolerance',1e-5,'FiniteDifferenceStepSize',0.1);
    
    % Reference = Driving
    if table.Reference(i) == 'Driving'
        drive_cost = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Distance(i);
              
        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * 5 + wait_time * 0 + drive_time * table.Driving_ref_1(i) + price * drive_cost + ASC_drive;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f1 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p1,1-p1,R) - CE;
        
        % Scenario 3
        CE = walk_time * 5 + wait_time * 0 + drive_time * table.Driving_ref_2(i) + price * drive_cost + ASC_drive; 
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f2 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p2,1-p2,R) - CE;
                      
        % Scenario 3
        CE = walk_time * 0 + wait_time * 0 + drive_time * table.Driving_ref_3(i) + price * drive_cost + ASC_drive;        
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 6 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        f3 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p3,1-p3,R) - CE;
        
        % Scenario 4
        CE = walk_time * 5 + wait_time * 0 + drive_time * table.Driving_ref_4(i) + price * drive_cost + ASC_drive;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 1 + wait_time * 1 + pool_time * 4 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        f4 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p4,1-p4,R) - CE;
        
        fun = @(x) system(f1,f2,f3,f4,x);
        cpt(i,:) = lsqnonlin(fun,[0.5,0.5,0.5,1.5],lb,ub,options);

%         cpt(i,:) = fsolve(fun,[0.5,0.5,0.5,1.5])
    elseif (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Subway (T)')
        
        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * transit_walk * 1.1 + wait_time * 9 + transit_time * table.Subway_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f1 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p1,1-p1,R) - CE;
        
        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 9 + transit_time * table.Subway_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f2 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p2,1-p2,R) - CE;
                      
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Subway_ref_3(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 6 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        f3 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p3,1-p3,R) - CE;
        
        % Scenario 4
        CE = walk_time * transit_walk + wait_time * 5 + transit_time * table.Subway_ref_4(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 1 + wait_time * 1 + pool_time * 4 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        f4 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p4,1-p4,R) - CE;
        
        fun = @(x) system(f1,f2,f3,f4,x);
        cpt(i,:) = lsqnonlin(fun,[0.5,0.5,0.5,1.5],lb,ub,options);

    elseif (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Bus')
        
        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * transit_walk * 1.1 + wait_time * 15 + transit_time * table.Bus_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f1 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p1,1-p1,R) - CE;
        
        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 15 + transit_time * table.Bus_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f2 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p2,1-p2,R) - CE;
                      
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Bus_ref_3(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 6 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        f3 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p3,1-p3,R) - CE;
        
        % Scenario 4
        CE = walk_time * transit_walk + wait_time * 8 + transit_time * table.Bus_ref_4(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 1 + wait_time * 1 + pool_time * 4 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        f4 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p4,1-p4,R) - CE;
        
        fun = @(x) system(f1,f2,f3,f4,x);
        cpt(i,:) = lsqnonlin(fun,[0.5,0.5,0.5,1.5],lb,ub,options);
        
    elseif (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Commuter rail')
        
        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * transit_walk * 1.1 + wait_time * 20 + transit_time * table.Rail_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f1 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p1,1-p1,R) - CE;
        
        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 20 + transit_time * table.Rail_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f2 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p2,1-p2,R) - CE;
                      
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 2 + transit_time * table.Rail_ref_3(i) + price * transit_cost + ASC_transit;        
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 6 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        f3 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p3,1-p3,R) - CE;
        
        % Scenario 4
        CE = walk_time * transit_walk + wait_time * 10 + transit_time * table.Rail_ref_4(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 1 + wait_time * 1 + pool_time * 4 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        f4 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p4,1-p4,R) - CE;
        
        fun = @(x) system(f1,f2,f3,f4,x);
        cpt(i,:) = lsqnonlin(fun,[0.5,0.5,0.5,1.5],lb,ub,options);
        
    elseif table.Reference(i) == 'Exclusive ridesharing'
        
        exc_cost = 2.2 + table.Exc_cost(i) * table.Distance(i) + 0.42 * table.Exc_driver_wait2(i);
                     
        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * 0 + wait_time * 9 + exc_time * 6 + price * table.Exc_ref_1(i) + ASC_exclusive;
        u1 = walk_time * 5 + wait_time * 2 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 4 + wait_time * 1 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f1 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p1,1-p1,R) - CE;
        
        % Scenario 2
        CE = walk_time * 0 + wait_time * 9 + exc_time * table.Exc_ref_2(i) + price * exc_cost * 1.2 + ASC_exclusive;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 5 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        f2 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p2,1-p2,R) - CE;
                      
        % Scenario 3
        CE = walk_time * 0 + wait_time * 1 + exc_time * table.Exc_ref_3(i) + price * exc_cost * 0.8 + ASC_exclusive;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 6 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        f3 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p3,1-p3,R) - CE;
        
        % Scenario 4
        CE = walk_time * 0 + wait_time * 5 + exc_time * table.Exc_ref_4(i) + price * exc_cost + ASC_exclusive;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 8 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 1 + wait_time * 1 + pool_time * 4 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        f4 = @(alpha,beta,gamma,lambda) subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p4,1-p4,R) - CE;
        
        fun = @(x) system(f1,f2,f3,f4,x);
        cpt(i,:) = lsqnonlin(fun,[0.5,0.5,0.5,1.5],lb,ub,options);
    end
end    

%% Population averages
alpha = mean(cpt(:,1))
beta = mean(cpt(:,2))
gamma = mean(cpt(:,3))
lambda = mean(cpt(:,4))

%% 
R = 0;
u = linspace(-100,100,10000);
p = linspace(0,1,10000);
V = zeros(1,length(u));
pi = zeros(1,length(p));

for i = 1:10000
    V(i) = value_fun(beta,gamma,lambda,R,u(i));
    pi(i) = exp(-(-log(p(i)))^alpha);
end

figure(1)
plot(u,V)
ylabel('Subjective value V');
xlabel('Objective utility');
title('Value function')

figure(2)
plot(p,pi)
ylabel('Subjective probability distortion');
xlabel('Objective probabilty');
title('Probability weighting function')

%%
function V = value_fun(beta,gamma,lambda,R,u)

if (u >= R)
    V = (u-R)^beta;
else
    V = -lambda*(R-u)^gamma;
end
end

%% 
function F = system(f1,f2,f3,f4,x)
% Here x(1) = alpha, x(2) = beta, x(3) = gamma, x(4) = lambda
F(1) = f1(x(1),x(2),x(3),x(4));
F(2) = f2(x(1),x(2),x(3),x(4));
F(3) = f3(x(1),x(2),x(3),x(4));
F(4) = f4(x(1),x(2),x(3),x(4));
end

%% Using Bernoulli CDF rather than PDF
function U = subjectiveU_v3(alpha,beta,gamma,lambda,u1,u2,p1,p2,R)
u = [u1,u2];
p_vals = [p1,p2];
[u_low,ind_low] = min(u);
[u_high,ind_high] = max(u);
p = p_vals(ind_low); 
if (u_high < R) % Pure loss
    w = exp(-(-log(p))^alpha);
    U = -w .* lambda .* (R - u_low).^gamma - (1-w) .* lambda .* (R - u_high).^gamma;
elseif (u_low >= R) % Pure gain
    w = exp(-(-log(1-p))^alpha);
    U = (u_low - R).^beta + w .* ((u_high - R).^beta - (u_low - R).^beta);
else % Mixed prospects: u_low < R & u_high > R
    w1 = exp(-(-log(p))^alpha); % gain
    w2 = exp(-(-log(1-p))^alpha); % loss
    U = -w1 .* lambda .* (R - u_low).^gamma + w2 .* (u_high - R).^beta;
end
end

%% Using new defn of weights
function U = subjectiveU_v2(alpha,beta,gamma,lambda,u1,u2,p1,p2,R)
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
function U = subjectiveU_v1(alpha,beta,gamma,lambda,u1,u2,R,pd)
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