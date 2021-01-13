%% PARAMETER ESTIMATION
% Vineet Jagadeesan Nair

% Using MATLAB's mnrfit function directly
% Assuming all coefficients to vary across modes

%% I - Generate data matrix

T = load('');
table_orig = T.Pilotdata117;

% Vary if table columns change
table = convertvars(table_orig,[2:5,19:26,32:116],'double');

num_responses = size(table,1); % Total no. of survey responses
scenarios = 5*num_responses; % 5 choice scenarios per respondent
Y = strings(scenarios,1); % Nominal response values
X = zeros(scenarios,13); % Observations for each scenario on 4 predictor variables

% Loop through each respondent
for i = 1:num_responses
    % Loop through each choice scenario
    for j = 1:5
        if table.transit_mode(i) == 'Bus'
            str1 = ["Mode_choice_bus_",num2str(j)]; 
            val1 = join(str1,"");
       elseif table.transit_mode(i) == 'Subway (T)'
           str1 = ["Mode_choice_subway_",num2str(j)]; 
           val1 = join(str1,"");
        end         
       
        Y(5*(i-1) + j) = table{i,val1};
            
        %% Public transit predictors
        % Walking times
        str2 = ["Transit_walk_variability_",num2str(j)]; 
        val2 = join(str2,"");
        X(5*(i-1) + j,1) = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * table{i,val2};

        if (table.transit_mode(i) == 'Subway (T)')
            % Waiting times
            str3 = ["Subway_wait_",num2str(j)]; 
            val3 = join(str3,"");
            X(5*(i-1) + j,2) = table{i,val3};
            
            % Riding times
            str4 = ["Subway_IVTT_",num2str(j)]; 
            val4 = join(str4,"");
            X(5*(i-1) + j,3) = table.Distance(i) * table{i,val4};
       
            % Price
            X(5*(i-1) + j,4) = table.Subway_cost(i);

        elseif (table.transit_mode(i) == 'Bus')
            % Waiting times
            str3 = ["Bus_wait_",num2str(j)]; 
            val3 = join(str3,"");
            X(5*(i-1) + j,2) = table{i,val3};
            
            % Riding times
            str4 = ["Bus_IVTT_",num2str(j)]; 
            val4 = join(str4,"");
            X(5*(i-1) + j,3) = table.Distance(i) * table{i,val4};
       
            % Price
            X(5*(i-1) + j,4) = table.Bus_cost(i)*rand;
        end
    
        %% Exclusive ridesharing predictors
        % Waiting times
        str5 = ["Exc_wait_",num2str(j)]; 
        val5 = join(str5,"");
        X(5*(i-1) + j,5) = table{i,val5};

        % Riding times
        str6 = ["Road_IVTT_",num2str(j)]; 
        val6 = join(str6,"");
        X(5*(i-1) + j,6) = table.Distance(i) * table{i,val6};
       
        % Price
        str7 = ["Exc_cost_variability",num2str(j)]; 
        val7 = join(str7,"");
        
        str8 = ["Exc_driver_wait",num2str(j)]; 
        val8 = join(str8,"");
        X(5*(i-1) + j,7) = table{i,val7} * (2.2 + 1.6 * table.Distance(i) + 0.42 * table{i,val8});
        
        %% Pooled ridesharing predictors
        % Walking times
        str9 = ["Pooled_walk_",num2str(j)]; 
        val9 = join(str9,"");
        X(5*(i-1) + j,8) = table{i,val9};

        % Waiting times
        str10 = ["Pooled_wait_",num2str(j)]; 
        val10 = join(str10,"");
        X(5*(i-1) + j,9) = table{i,val10};

        % Riding times
        str11 = ["Pooled_IVTT",num2str(j)]; 
        val11 = join(str11,"");
        X(5*(i-1) + j,10) = table.Distance(i) + table{i,val11};
       
        % Price
        X(5*(i-1) + j,11) = table{i,val7} * (2.2 + table.Pool_cost(i) * table.Distance(i));
        
        %% Driving predictors - reference/baseline alternative
        % Walking times
        str12 = ["Road_walk_",num2str(j)]; 
        val12 = join(str12,"");
        X(5*(i-1) + j,12) = table{i,val12};
        
        % Riding tie with personal car same as exclusive rideshare
        % So don't include separately - would be linearly dependent column!
       
        % Price
        str13 = ["Road_cost_",num2str(j)]; 
        val13 = join(str13,"");
        X(5*(i-1) + j,13) = table.Parking_cost(i) + table.Cost_per_mile(i) * table{i,val13} * table.Distance(i);
    end
end

%% II - Estimation of mode-choice model
% Multinomial logit for utility functions

Y_resp = categorical(Y);
[B,dev,stats] = mnrfit(X,Y_resp);

% Standard errors
stats.se

% Statistical significance (p-values)
stats.p

% 95% confidence intervals for coefficients
LL = stats.beta - 1.96.*stats.se;
UL = stats.beta + 1.96.*stats.se;

% Deviance of fit 
% Difference between maximum achievable log likelihood and that attained under fitted model
dev
