%% CPT parameters derivation

% Using originally intended approach of estimating using parameters separately using pure gain, pure loss and mixed

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
% Need to adjust regression to account for wages 
% Annual income --> Hourly wage

VOT_walking = abs(price/walk_time)*100
VOT_waiting = abs(price/wait_time)*100
VOT_drive = abs(price/drive_time)*100
VOT_transit = abs(price/transit_time)*100
VOT_rideshare = abs(price/rideshare_time)*100

%% Risk scenarios

T = load('Pilotdata1218');
table_orig = T.Pilotdata1218;
% May need to modify column indices later on
table = convertvars(table_orig,[2:5,14:21,27:111],'double');

%%

num_responses = size(table,1); % Total no. of survey responses
scenarios = 5*num_responses; % 5 choice scenarios per respondent
alternatives = scenarios*4; % 4 alternatives per choice situation
X = zeros(alternatives,13); % Observations for each scenario on 4 predictor variables

% Each row gives risk attitude parameters - alpha, beta+, beta-, lambda
cpt = zeros(num_responses,4);

for i = 8%1:num_responses
    
    % CDF of SMoDS objective utility - assume extreme value distribution (Gumbel - type I)
    % Could also use truncated normal
    % Best and worst case utilities - correspond to shortest & longest travel times
    price_low = (2.2 + table.Pool_cost(i)*table.Distance(i))*0.8;
    price_high = (2.2 + table.Pool_cost(i)*table.Distance(i))*1.2;
    u_best = walk_time * 0 + wait_time * 1 + rideshare_time * 6 * table.Distance(i) + price * price_low + ASC_pooled;
    u_worst = walk_time * 0 + wait_time * 1 + rideshare_time * 6 * table.Distance(i) + price * price_high + ASC_pooled;
    
    mu = (u_best + u_worst)/2; % Location parameter
    sigma = u_best - u_worst; % Scale parameter
    pd = makedist('Extreme Value','mu',mu,'sigma',sigma);
    
    % Reference = Driving
    if table.Reference(i) == 'Driving'
        
        ref_cost = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Distance(i);

        % Scenario 1: Pure-gain 1
        R = walk_time * 5 + wait_time * 0 + drive_time * table.Driving_ref_1(i) + price * ref_cost + ASC_drive;
        
        SMODS_cost = (2.2 + table.Pool_cost(i)*table.Distance(i))*0.8; 
        u1 = walk_time * 2 + wait_time * 0 + rideshare_time * 7 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.1;

        u2 = walk_time * 3 + wait_time * 1 + rideshare_time * 6 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.9;
        
        f1 = @(alpha,beta) puregain(alpha,beta,u1,u2,R,pd);
        
        % Scenario 2: Pure gain 2
        R = walk_time * 5 + wait_time * 0 + drive_time * table.Driving_ref_2(i) + price * ref_cost + ASC_drive;
        p1 = 0.6;
        p2 = 0.4;
        
        f2 = @(alpha,beta) puregain(alpha,beta,u1,u2,R,pd);
        fun1 = @(x) system(f1,f2,x);
%         x = fsolve(fun1,[0.5,0.5])
        
        % Non-linear least squares curve fitting
        lb = [0.01,0.01];
        ub = [0.99,0.99];
        options = optimoptions(@lsqnonlin,'Display','iter','OptimalityTolerance',1e-10, ...
            'DiffMinChange',0.1,'FiniteDifferenceStepSize',1e-1);
        cpt(i,1:2) = lsqnonlin(fun1,[0.5,0.5],lb,ub,options)

        % Trying to find global minimum
%         f_val = @(x) fun1(x);
%         squared_error = @(x) f_val(1).*f_val(1) + f_val(2).*f_val(2);
%         opts = optimoptions(@fmincon,'Algorithm','active-set');
%         problem = createOptimProblem('fmincon','x0',[0.5,0.5],'objective',squared_error, ...
%             'lb',[1e-6;1e-6],'ub',[1;1],'options',opts);
%         gs = GlobalSearch;
%         rng(14,'twister') % for reproducibility
%         [xfinal,fval] = run(gs,problem)

        % Scenario 3: Pure loss
        SMODS_cost = (2.2 + table.Pool_cost(i)*table.Distance(i))*1.2; 
        
        R = walk_time * 0 + wait_time * 0 + drive_time * table.Driving_ref_3(i) + price * ref_cost + ASC_drive;
        u1 = walk_time * 7 + wait_time * 4 + rideshare_time * 6 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.8;

        u2 = walk_time * 5 + wait_time * 6 + rideshare_time * 10 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.2;
        
        alpha = cpt(i,1);
        beta = cpt(i,2);
        f3 = @(gamma,lambda) pureloss(alpha,gamma,lambda,u1,u2,R,pd);
        
        % Scenario 4: Mixed prospects
        SMODS_cost = 2.2 + table.Pool_cost(i)*table.Distance(i); 
        
        R = walk_time * 0 + wait_time * 0 + drive_time * table.Driving_ref_3(i) + price * ref_cost + ASC_drive;
        u1 = walk_time * 7 + wait_time * 6 + rideshare_time * 9 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p1 = 0.5; % loss

        u2 = walk_time * 1 + wait_time * 1 + rideshare_time * 7 * table.Distance(i) + price * SMODS_cost + ASC_pooled;
        p2 = 0.5; % gain
        
        alpha = cpt(i,1);
        f4 = @(gamma,lambda) mixed(alpha,beta,gamma,lambda,u1,u2,R,pd);
        fun2 = @(y) system(f3,f4,y);
        
        lb = [0.01,1];
        ub = [0.99,inf];
        options = optimoptions(@lsqnonlin,'Display','iter','OptimalityTolerance',1e-10, ...
            'DiffMinChange',0.1,'FiniteDifferenceStepSize',1e-1);
        cpt(i,3:4) = lsqnonlin(fun2,[0.5,0.5],lb,ub,options)
    end
end        
        
%% Plotting function
alpha = 0:0.1:1;
beta = 0:0.1:1;
error = zeros(length(alpha),length(beta));
for i = 1:length(alpha)
    for j = 1:length(beta)
        x = [alpha(i),beta(j)];
        error(i,j) = norm(fun(x));
    end
end
contourf(alpha,beta,error)

%% 
function F = system(f1,f2,x)
% Here x(1) = alpha & x(2) = beta OR x(1) = gamma & x(2) = lambda
F(1) = f1(x(1),x(2));
F(2) = f2(x(1),x(2));
end

% Where do p1 and p2 even play a role here??
function F = puregain(alpha,beta,u1,u2,R,pd)
w1 = -exp(-(-log(1-cdf(pd,u1)))^alpha);
w2 = exp(-(-log(1-cdf(pd,u1)))^alpha) - exp(-(-log(1-cdf(pd,u2)))^alpha);
F = w1 .* abs(u1 - R).^beta + w2 .* abs(u2 - R).^beta - R;
end

function F = pureloss(alpha,gamma,lambda,u1,u2,R,pd)
w1 = exp(-(-log(cdf(pd,u1)))^alpha);
w2 = exp(-(-log(cdf(pd,u2)))^alpha) - exp(-(-log(cdf(pd,u1)))^alpha);
F = -w1 .* lambda .* abs(R - u1).^gamma - w2 .* lambda .* abs(R - u2).^gamma - R;
end

function F = mixed(alpha,beta,gamma,lambda,u1,u2,R,pd)
w1 = exp(-(-log(cdf(pd,u1)))^alpha); % loss
w2 = exp(-(-log(1-cdf(pd,u1)))^alpha) - exp(-(-log(1-cdf(pd,u2)))^alpha); % gain
F = -w1 .* lambda .* abs(R - u1).^gamma + w2 .* abs(u2 - R).^beta - R;
end

    