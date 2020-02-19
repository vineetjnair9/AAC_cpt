%% PARAMETER ESTIMATION
% Vineet Jagadeesan Nair
% 19 February 2020

% Maximum simulated likelihood estimation by Kenneth Train

% Updated for new study design
% Assumptions:
% (1) All coefficients as fixed not random
% 

%% I - Estimation of mode-choice model
% Multinomial logit for utility functions
% Generation of data matrix X for using Kenneth Train's code

T = load('Pilotdata0117');
table_orig = T.Pilotdata117;

% Vary if table columns change
table = convertvars(table_orig,[2:5,19:26,32:116],'double');

num_responses = size(table,1); % Total no. of survey responses
scenarios = 5*num_responses; % 5 choice scenarios per respondent
alternatives = scenarios*4; % 4 alternatives per choice situation

% Vary if table columns change
X = zeros(alternatives,13); % Observations for each scenario on 4 predictor variables

row = 1;
choice = 1;

% Loop through each respondent
for i = 1:num_responses
    % Loop through each choice scenario
    for j = 1:5
        % Loop through each alternative
        for k = 1:4
            X(row,1) = i; % Person facing this alternative
            X(row,2) = choice; % Choice situations numbered sequentially
        
           if table.transit_mode(i) == 'Bus'
                str1 = ["Mode_choice_bus_",num2str(j)]; 
                val1 = join(str1,"");
           elseif table.transit_mode(i) == 'Subway (T)'
               str1 = ["Mode_choice_subway_",num2str(j)]; 
               val1 = join(str1,"");
           end                
            
            if k == 1 % Driving alternative
                if(table{i,val1} == 'Driving')
                    X(row,3) = 1;
                end
                            
                % Walking times
                str2 = ["Road_walk_",num2str(j)]; 
                val2 = join(str2,"");
                X(row,4) = table{i,val2};

                % Waiting times 
                % Either 0 min (but this makes X singular) or drawn from a uniform random dist between 0-1 min
                X(row,5) = 0;

                % Riding times
                str3 = ["Road_IVTT_",num2str(j)]; 
                val3 = join(str3,"");
                X(row,6) = table.Distance(i) * table{i,val3};

                % Price
                str4 = ["Road_cost_",num2str(j)]; 
                val4 = join(str4,"");
                X(row,9) = table.Parking_cost(i) + table.Cost_per_mile(i) * table{i,val4} * table.Distance(i);

                X(row,10) = 1;
        
            elseif k == 2 % Public transit alternative
                if(table{i,val1} == 'Public transit')
                    X(row,3) = 1;
                end    
                
                % Walking times
                str5 = ["Transit_walk_variability_",num2str(j)]; 
                val5 = join(str5,"");
                X(row,4) = (table.Walk_hometransit(i) + table.Walk_transitwork(i)) * table{i,val5};

                if (table.transit_mode(i) == 'Subway (T)')
                    % Waiting times
                    str6 = ["Subway_wait_",num2str(j)]; 
                    val6 = join(str6,"");
                    X(row,5) = table{i,val6};

                    % Riding times
                    str7 = ["Subway_IVTT_",num2str(j)]; 
                    val7 = join(str7,"");
                    X(row,7) = table.Distance(i) * table{i,val7};

                    % Price
                    X(row,9) = table.Subway_cost(i);

                elseif (table.transit_mode(i) == 'Bus')
                    % Waiting times
                    str8 = ["Bus_wait_",num2str(j)]; 
                    val8 = join(str8,"");
                    X(row,5) = table{i,val8};

                    % Riding times
                    str9 = ["Bus_IVTT_",num2str(j)]; 
                    val9 = join(str9,"");
                    X(row,7) = table.Distance(i) * table{i,val9};

                    % Price
                    X(row,9) = table.Bus_cost(i);
                end
                
                X(row,11) = 1;
                
            elseif k == 3 % Exclusive rideshare alternative
                if(table{i,val1} == 'Exclusive rideshare')
                    X(row,3) = 1;
                end   
                
                % Walking times
                % Either 0 min (but this makes X singular) or drawn from a uniform random dist between 0-1 min
                X(row,4) = 0;

                % Waiting times
                str10 = ["Exc_wait_",num2str(j)]; 
                val10 = join(str10,"");
                X(row,5) = table{i,val10};

                % Riding times
                X(row,8) = table.Distance(i) * table{i,val3};

                % Price
                str11 = ["Exc_cost_variability",num2str(j)]; 
                val11 = join(str11,"");

                str12 = ["Exc_driver_wait",num2str(j)]; 
                val12 = join(str12,"");
                X(row,9) = table{i,val11} * (2.2 + 1.6 * table.Distance(i) + 0.42 * table{i,val12});
        
                X(row,12) = 1;
                
            elseif k == 4 % Pooled rideshare alternative
                if(table{i,val1} == 'Pooled rideshare')
                    X(row,3) = 1;
                end   
                
                % Walking times
                str13 = ["Pooled_walk_",num2str(j)]; 
                val13 = join(str13,"");
                X(row,4) = table{i,val13};

                % Waiting times
                str14 = ["Pooled_wait_",num2str(j)]; 
                val14 = join(str14,"");
                X(row,5) = table{i,val14};

                % Riding times
                str15 = ["Pooled_IVTT",num2str(j)]; 
                val15 = join(str15,"");
                X(row,8) = table.Distance(i) * table{i,val15};

                % Price
                X(row,9) = table{i,val11} * (2.2 + table.Pool_cost(i) * table.Distance(i));
                X(row,13) = 1;
            end
            row = row + 1;
        end 
        choice = choice + 1;
    end
end

writematrix(X,'data_new.txt');
    