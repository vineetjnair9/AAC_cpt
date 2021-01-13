%% CPT parameters derivation
clear; clc;

% Mode-choice parameters
load('paramhat.mat')
load('valid_indices.mat')

% Declare GLOBAL variables
global walk_time wait_time transit_time exc_time pool_time
global price ASC_transit ASC_exclusive ASC_pooled

walk_time = -paramhat(1);
wait_time = -paramhat(2);
transit_time = -paramhat(3);
exc_time = -paramhat(4);
pool_time = -paramhat(5);
price = -paramhat(6);
ASC_transit = 0; % Baseline
ASC_exclusive = paramhat(7);
ASC_pooled = paramhat(8);

disp('Value of time in $/h:')
VOT_walking = (walk_time/price)*60
VOT_waiting = (wait_time/price)*60
VOT_transit = (transit_time/price)*60
VOT_exc = (exc_time/price)*60
VOT_pool = (pool_time/price)*60

% Risk scenarios

T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_6-10_n955');
table_orig = T.Fulllaunch610n955;

% Vary if table columns change
table = convertvars(table_orig,[7:8,10:12,48:139],'double');
num_responses = size(table,1); % Total no. of survey responses

% Each row gives risk attitude parameters - alpha, beta+, beta-, lambda
cpt = zeros(num_responses,5);
error = zeros(num_responses,1); % Squared 2-norm of residual

panel = @(x,respondent_num) panel_obj(x,respondent_num,table,'static_SMoDS','modified',paramhat);

% Settings for non linear least squares
lb = [0,0,0,0,0];
ub = [1,1,1,1,1];
x0 = [0.5,0.5,0.5,0.5,0.5];

options = optimoptions('lsqnonlin','Display','off','FiniteDifferenceStepSize',1e-2,'DiffMinChange',1e-2,'DiffMaxChange',0.1,...
    'FiniteDifferenceType','central','MaxFunctionEvaluations',5000,'MaxIterations',4000);
%'UseParallel',true,'FunctionTolerance',1e-10,'OptimalityTolerance',1e-10

%% Estimating CPT parameters on all responses
tic
for i = 1:20
    fun = @(x) panel(x,i);
    [cpt(i,:),error(i)] = lsqnonlin(fun,x0,lb,ub,options);
end
toc

%% Estimating CPT parameters only using valid responses
tic
for i = 1:20%length(valid_indices)   
    j = valid_indices(i);
    fun = @(x) panel(x,j);
    [cpt(i,:),error(i)] = lsqnonlin(fun,x0,lb,ub,options);
end
toc

% R_type = dynamic_u0, dynamic_SMoDS, dynamic_mean_prob, dynamic_mean, static_u0, static_SMoDS
% weighting = original, modified
% optimization method = 'lsq' (nonlinear least squares), 'grid' (grid search), 'newton', 'global-opt' (Global optimization toolbox)

% With normalized lambda
cpt(:,5) = cpt(:,5) * 100;

%% Objective function for each combination of parameter values under consideration
function F = panel_obj(x,respondent_num,table,R_type,weight_type,paramhat)
    
    % Normalize lambda so that all 4 parameters lie between 0 and 1
    x(5) = x(5) * 100;

    walk_time = -paramhat(1);
    wait_time = -paramhat(2);
    transit_time = -paramhat(3);
    exc_time = -paramhat(4);
    pool_time = -paramhat(5);
    price = -paramhat(6);
    ASC_transit = 0; % Baseline
    ASC_exclusive = paramhat(7);
    ASC_pooled = paramhat(8);
       
    ref = @(u1,u2,p,CE,i) reference(u1,u2,p,R_type,CE,table,i,paramhat);
    calc_error = @(a,b,c,d,e,CE_sR,u1,u2,p,R) error_func(a,b,c,d,e,CE_sR,u1,u2,p,R,weight_type);
    
    i = respondent_num;
    transit_walk = table.transit_walk_1(i) + table.transit_walk_2(i);
    transit_cost = table.Transit_cost(i);

    SMODS_cost = 2.2 + table.Pool_cost(i)*table.Distance_corrected(i);

    p1 = table.p1(i)/100;
    p2 = table.p2(i)/100;
    p3 = table.p3(i)/100;
    p4 = table.p4(i)/100;
    p5 = table.p5(i)/100;
    p6 = table.p6(i)/100;

    if (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Bus')

        % Scenario 1
        CE = walk_time * transit_walk * 1.1 + wait_time * 15 + transit_time * table.Bus_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p1,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p1,R);

        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 15 + transit_time * table.Bus_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p2,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p2,R);
        
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Bus_ref_3(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 6.5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p3,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p3,R);
        
        % Scenario 4     
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Bus_ref_4(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 6 + wait_time * 5 + pool_time * 7 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 9 + wait_time * 4 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p4,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p4,R);
        
        % Scenario 5
        CE = walk_time * transit_walk + wait_time * 8 + transit_time * table.Bus_ref_5(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 2 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p5,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p5,R);
        
        % Scenario 6
        CE = walk_time * transit_walk + wait_time * 5 + transit_time * table.Bus_ref_6(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 5 + wait_time * 6 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p6,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p6,R);
        
        F = [f1; f2; f3; f4; f5; f6]; 

    elseif (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Subway (T)')

        % Scenario 1
        CE = walk_time * transit_walk * 1.1 + wait_time * 9 + transit_time * table.Subway_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p1,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p1,R);
        
        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 9 + transit_time * table.Subway_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 2 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 1 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p2,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p2,R);
        
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Subway_ref_3(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 8 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p3,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p3,R);
        
        % Scenario 4  
        CE = walk_time * transit_walk * 0.9 + wait_time * 1 + transit_time * table.Subway_ref_4(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 6 + wait_time * 5 + pool_time * 7 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 9 + wait_time * 4 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p4,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p4,R);
        
        % Scenario 5
        CE = walk_time * transit_walk + wait_time * 4 + transit_time * table.Subway_ref_5(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 2 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p5,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p5,R);
        
        % Scenario 6
        CE = walk_time * transit_walk + wait_time * 5 + transit_time * table.Subway_ref_6(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 5 + wait_time * 6 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p6,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p6,R);
        
        F = [f1; f2; f3; f4; f5; f6];

    elseif (table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}' && table.transit_mode(i) == 'Commuter rail')

        % Scenario 1
        CE = walk_time * transit_walk * 1.1 + wait_time * 15 + transit_time * table.Rail_ref_1(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 4 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 3 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p1,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p1,R);
        
        % Scenario 2
        CE = walk_time * transit_walk * 1.1 + wait_time * 20 + transit_time * table.Rail_ref_2(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 2 + wait_time * 4 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p2,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p2,R);
        
        % Scenario 3      
        CE = walk_time * transit_walk * 0.9 + wait_time * 2 + transit_time * table.Rail_ref_3(i) + price * transit_cost + ASC_transit;        
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 6.5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p3,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p3,R);
        
        % Scenario 4      
        CE = walk_time * transit_walk * 0.9 + wait_time * 2 + transit_time * table.Rail_ref_4(i) + price * transit_cost + ASC_transit;        
        u1 = walk_time * 6 + wait_time * 5 + pool_time * 7 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 9 + wait_time * 4 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p4,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p4,R);
        
        % Scenario 5
        CE = walk_time * transit_walk + wait_time * 10 + transit_time * table.Rail_ref_5(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 2 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p5,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p5,R);
        
        % Scenario 6
        CE = walk_time * transit_walk + wait_time * 10 + transit_time * table.Rail_ref_6(i) + price * transit_cost + ASC_transit;
        u1 = walk_time * 5 + wait_time * 6 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p6,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p6,R);
        
        F = [f1; f2; f3; f4; f5; f6];

    elseif table.Reference(i) == 'Exclusive ridesharing'

        exc_cost = 2.2 + table.Exc_cost(i) * table.Distance_corrected(i) + 0.42 * table.Exc_driver_wait2(i);

        % Scenario 1
        % Certainty equivalent - Sure prospect with certain alternative travel options
        CE = walk_time * 0 + wait_time * 9 + exc_time * table.Distance_corrected(i) * 5 + price * table.Exc_ref_1(i) + ASC_exclusive;
        u1 = walk_time * 5 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 4 + wait_time * 1 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p1,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f1 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p1,R);
        
        % Scenario 2
        CE = walk_time * 0 + wait_time * 9 + exc_time * table.Exc_ref_2(i) + price * exc_cost * 1.2 + ASC_exclusive;
        u1 = walk_time * 2 + wait_time * 4 + pool_time * 4.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 3 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost * 0.8 + ASC_pooled;
        R = ref(u1,u2,p2,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f2 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p2,R);
        
        % Scenario 3
        CE = walk_time * 0 + wait_time * 1 + exc_time * table.Exc_ref_3(i) + price * exc_cost * 0.8 + ASC_exclusive;
        u1 = walk_time * 7 + wait_time * 4 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 5 + wait_time * 6 + pool_time * 6.5 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p3,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f3 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p3,R);
        
        % Scenario 4
        CE = walk_time * 0 + wait_time * 1 + exc_time * table.Exc_ref_4(i) + price * exc_cost * 0.8 + ASC_exclusive;
        u1 = walk_time * 6 + wait_time * 5 + pool_time * 7 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        u2 = walk_time * 9 + wait_time * 4 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost * 1.2 + ASC_pooled;
        R = ref(u1,u2,p4,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f4 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p4,R);
        
        % Scenario 5
        CE = walk_time * 0 + wait_time * 5 + exc_time * table.Exc_ref_5(i) + price * exc_cost + ASC_exclusive;
        u1 = walk_time * 7 + wait_time * 6 + pool_time * 6 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 2 + wait_time * 3 + pool_time * 3 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p5,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f5 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p5,R);
        
        % Scenario 6
        CE = walk_time * 0 + wait_time * 5 + exc_time * table.Exc_ref_6(i) + price * exc_cost + ASC_exclusive;
        u1 = walk_time * 5 + wait_time * 6 + pool_time * 5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        u2 = walk_time * 3 + wait_time * 2 + pool_time * 3.5 * table.Distance_corrected(i) + price * SMODS_cost + ASC_pooled;
        R = ref(u1,u2,p6,CE,i);
%         CE_sR = CE; 
        CE_sR = value_func(x(3),x(4),x(5),CE,R);
        f6 = calc_error(x(1),x(2),x(3),x(4),x(5),CE_sR,u1,u2,p6,R);
        
        F = [f1; f2; f3; f4; f5; f6];
       
    end
end

%% Calculate error for nonlinear equation
function error = error_func(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,CE_actual,u1,u2,p,R,weight_type)
    CE_pred = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p,R,weight_type);
    
    % No normalization 
%     error = CE_pred - CE_actual;
    
    % With normalization
    if (CE_actual == 0)
        error = abs(CE_pred - CE_actual);
    else         
        % Taking relative error w.r.t true value (CE_sR)
        error = (CE_pred - CE_actual)/CE_actual;

      % Taking relative error w.r.t max subjective value
%         error = (CE_pred - CE_actual)/max([abs(CE_pred),abs(CE_actual)]);        
    end

end

%% Calculating subjective utilities
function U_sR = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p1,R,weight_type) 

    V = @(u) value_func(beta_plus,beta_minus,lambda,u,R);
    u = [u1,u2];
    p_vals = [p1,1-p1];

    [u_low,ind_low] = min(u);
    [u_high,~] = max(u);

    p = p_vals(ind_low);

    w = weights(p,u_low,u_high,R,alpha_plus,alpha_minus,weight_type);
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
        
%% Probability weighting function
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
%% Determining reference point
function R = reference(u1,u2,p,R_type,CE,table,i,paramhat)
    walk_time = -paramhat(1);
    wait_time = -paramhat(2);
    transit_time = -paramhat(3);
    exc_time = -paramhat(4);
    pool_time = -paramhat(5);
    price = -paramhat(6);
    ASC_transit = 0; % Baseline
    ASC_exclusive = paramhat(7);
    ASC_pooled = paramhat(8);
    
    if strcmp(R_type,'dynamic_u0') % Setting R as the certainty equivalent itself 
        % i.e. R = objective utility of current most frequent alternative
        R = CE;
    elseif strcmp(R_type,'dynamic_SMoDS')
        % R = mean of 2 SMoDS outcomes
        R = p*u1 + (1-p)*u2;
    elseif strcmp(R_type,'dynamic_mean_prob')
        % R = average of current alternative and mean SMoDS outcome
        % Accounting for probabilities
        R = (CE + p*u1 + (1-p)*u2)/2;  
    elseif strcmp(R_type,'dynamic_mean')
        % R = simple mean of all 3 outcomes (no probabilities)
        R = (CE + u1 + u2)/3;
    elseif strcmp(R_type,'static_u0')
        if table.Reference(i) == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}'
            R = walk_time * (table.transit_walk_1(i) + table.transit_walk_2(i)) + wait_time * 5 + transit_time * 3 * table.Distance(i) + price * table.Transit_cost(i) + ASC_transit;
        else
            Exc_cost = 2.2 + table.Exc_cost(i)*table.Distance(i);
            R = walk_time * 0 + wait_time * 3.5 + exc_time * 4 * table.Distance(i) + price * Exc_cost + ASC_exclusive;
        end
    elseif strcmp(R_type,'static_SMoDS')
        SMoDS_cost_base = 2.2 + table.Pool_cost(i)*table.Distance(i);
        R = walk_time * 3.5 + wait_time * 3.5 + pool_time * 5 * table.Distance(i) + price * SMoDS_cost_base + ASC_pooled;
    end
 
end