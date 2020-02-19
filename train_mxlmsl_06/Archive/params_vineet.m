%% PARAMETER ESTIMATION
% Vineet Jagadeesan Nair
% 17 December 2019

%% I - Estimation of mode-choice model
% Multinomial logit for utility functions

T = load('Pilot_data_12-18');
table = T.data;

num_responses = size(table,1); % Total no. of survey responses
scenarios = 5*num_responses; % 5 choice scenarios per respondent
Y = zeros(scenarios,1); % Nominal response values
X = zeros(scenarios,4); % Observations for each scenario on 4 predictor variables

% Loop through each respondent
for i = 1:num_responses
    Y(5*(i-1) + 1) = table.Mode_choice_1(i);
    Y(5*(i-1) + 2) = table.Mode_choice_2(i);
    Y(5*(i-1) + 3) = table.Mode_choice_3(i);
    Y(5*(i-1) + 4) = table.Mode_choice_4(i);
    Y(5*(i-1) + 5) = table.Mode_choice_5(i);
    
    %% Driving predictors
    % Walking times
    X(5*(i-1) + 1,1) = table.Road_walk_1(i);
    X(5*(i-1) + 2,1) = table.Road_walk_2(i);
    X(5*(i-1) + 3,1) = table.Road_walk_3(i);
    X(5*(i-1) + 4,1) = table.Road_walk_4(i);
    X(5*(i-1) + 5,1) = table.Road_walk_5(i);
    
    % Waiting times
    X(5*(i-1) + 1,2) = 0;
    X(5*(i-1) + 2,2) = 0;
    X(5*(i-1) + 3,2) = 0;
    X(5*(i-1) + 4,2) = 0;
    X(5*(i-1) + 5,2) = 0;
    
    % Riding times
    X(5*(i-1) + 1,3) = table.Distance(i) * table.Road_IVTT_1(i);
    X(5*(i-1) + 2,3) = table.Distance(i) * table.Road_IVTT_2(i);
    X(5*(i-1) + 3,3) = table.Distance(i) * table.Road_IVTT_3(i);
    X(5*(i-1) + 4,3) = table.Distance(i) * table.Road_IVTT_4(i);
    X(5*(i-1) + 5,3) = table.Distance(i) * table.Road_IVTT_5(i);
    
    % Price
    X(5*(i-1) + 1,4) = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Road_cost_1(i) * table.Distance(i);
    X(5*(i-1) + 2,4) = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Road_cost_2(i) * table.Distance(i);
    X(5*(i-1) + 3,4) = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Road_cost_3(i) * table.Distance(i);
    X(5*(i-1) + 4,4) = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Road_cost_4(i) * table.Distance(i);
    X(5*(i-1) + 5,4) = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Road_cost_5(i) * table.Distance(i);
    
    %% Public transit predictors
    % Walking times
    X(5*(i-1) + 1,5) = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * table.Transit_walk_variability_1(i);
    X(5*(i-1) + 2,5) = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * table.Transit_walk_variability_2(i);
    X(5*(i-1) + 3,5) = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * table.Transit_walk_variability_3(i);
    X(5*(i-1) + 4,5) = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * table.Transit_walk_variability_4(i);
    X(5*(i-1) + 5,5) = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * table.Transit_walk_variability_5(i);
    
    if (table.transit_mode(i) == 'Subway (T)')
        % Waiting times
        X(5*(i-1) + 1,6) = table.Subway_wait_1(i);
        X(5*(i-1) + 2,6) = table.Subway_wait_2(i);
        X(5*(i-1) + 3,6) = table.Subway_wait_3(i);
        X(5*(i-1) + 4,6) = table.Subway_wait_4(i);
        X(5*(i-1) + 5,6) = table.Subway_wait_5(i);

        % Riding times
        X(5*(i-1) + 1,7) = table.Distance(i) * table.Subway_IVTT_1(i);
        X(5*(i-1) + 2,7) = table.Distance(i) * table.Subway_IVTT_2(i);
        X(5*(i-1) + 3,7) = table.Distance(i) * table.Subway_IVTT_3(i);
        X(5*(i-1) + 4,7) = table.Distance(i) * table.Subway_IVTT_4(i);
        X(5*(i-1) + 5,7) = table.Distance(i) * table.Subway_IVTT_5(i);

        % Price
       for j = 1:5
            X(5*(i-1) + j,8) = table.Subway_cost(i);
        end
        
    elseif (table.transit_mode(i) == 'Bus')
        % Waiting times
        X(5*(i-1) + 1,6) = table.Bus_wait_1(i);
        X(5*(i-1) + 2,6) = table.Bus_wait_2(i);
        X(5*(i-1) + 3,6) = table.Bus_wait_3(i);
        X(5*(i-1) + 4,6) = table.Bus_wait_4(i);
        X(5*(i-1) + 5,6) = table.Bus_wait_5(i);

        % Riding times
        X(5*(i-1) + 1,7) = table.Distance(i) * table.Bus_IVTT_1(i);
        X(5*(i-1) + 2,7) = table.Distance(i) * table.Bus_IVTT_2(i);
        X(5*(i-1) + 3,7) = table.Distance(i) * table.Bus_IVTT_3(i);
        X(5*(i-1) + 4,7) = table.Distance(i) * table.Bus_IVTT_4(i);
        X(5*(i-1) + 5,7) = table.Distance(i) * table.Bus_IVTT_5(i);

        % Price
        for j = 1:5
            X(5*(i-1) + j,8) = table.Bus_cost(i);
        end
    end
    
    %% Exclusive ridesharing predictors
    