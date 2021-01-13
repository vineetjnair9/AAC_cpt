%% PARAMETER ESTIMATION
% Vineet Jagadeesan Nair
% 19 February 2020

% Maximum simulated likelihood estimation by Kenneth Train

% Updated for new study design
% Assumptions:
% (1) All coefficients as fixed not random
% (2) Using Bernoulli CDF

%% I - Estimation of mode-choice model
% Multinomial logit for utility functions
% Generation of data matrix X for using Kenneth Train's code

%data = xlsread('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Pilot_v2_2-20-20');

T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Pilot_v2_2-20-20.mat');

load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Code\AAC_cpt\Survey analysis\factor_levels.mat')

table_orig = T.Pilotv222020;

% Vary if table columns change
table = convertvars(table_orig,[1,4:5,7:8,43:125],'double');

num_responses = size(table,1); % Total no. of survey responses
scenarios = 5*num_responses; % 5 choice scenarios per respondent
alternatives = scenarios*4; % 4 alternatives per choice situation

X = zeros(alternatives,14); % Observations for each scenario on 4 predictor variables            
row = 1;
choice = 1;

% Loop through each respondent
for i = 1:num_responses
    
    % Finding which scenarios are randomly drawn for each respondent
    choice_sets = find(~cellfun(@isempty,table2array(table_orig(i,10:42))));
    
    if table.transit_mode(i) == 'Subway (T)'
            choice_sets = choice_sets - 11;
    elseif table.transit_mode(i) == 'Commuter rail'
            choice_sets = choice_sets - 22;
    end    
    
    % Loop through each choice scenario
    for j = 1:5
        
        if table.transit_mode(i) == 'Bus'    
            str1 = ["bus",num2str(choice_sets(j))]; 
            val1 = join(str1,"");
        elseif table.transit_mode(i) == 'Subway (T)'
            str1 = ["subway",num2str(choice_sets(j))]; 
            val1 = join(str1,"");
        elseif table.transit_mode(i) == 'Commuter rail'
            str1 = ["Rail",num2str(choice_sets(j))]; 
            val1 = join(str1,"");
        end    
        
        % Find correct factor levels
        walk = factor_levels(choice_sets(j),1);
        wait = factor_levels(choice_sets(j),2);
        exc_driver_wait = factor_levels(choice_sets(j),3);
        drive_ride = factor_levels(choice_sets(j),4);
        transit_ride = factor_levels(choice_sets(j),5);
        exc_ride = factor_levels(choice_sets(j),6);
        pool_ride = factor_levels(choice_sets(j),7);
        exc_cost_var = factor_levels(choice_sets(j),8);
        pool_cost_var = factor_levels(choice_sets(j),9);
        
        % Loop through each alternative
        for k = 1:4
            X(row,1) = i; % Person facing this alternative
            X(row,2) = choice; % Choice situations numbered sequentially            
            
            if k == 1 % Driving alternative
                if(table{i,val1} == 'Driving')
                    X(row,3) = 1;
                end
                            
                % Walking times
                str = ["Road_walk_",num2str(walk)]; 
                val = join(str,"");
                X(row,4) = table{i,val};

                % Waiting times 
                % Either 0 min (but this makes X singular) or drawn from a uniform random dist between 0-1 min
                X(row,5) = 0;

                % Riding times
                str = ["Road_IVTT_",num2str(drive_ride)]; 
                val = join(str,"");
                X(row,6) = table.Distance(i) * table{i,val};

                % Price
                X(row,10) = table.Parking_cost(i) + table.Cost_per_mile(i) * table.Distance(i);

                X(row,11) = 1;
        
            elseif k == 2 % Public transit alternative
                if(table{i,val1} == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}')
                    X(row,3) = 1;
                end    
                
                % Walking times
                str = ["Transit_walk_variability_",num2str(walk)]; 
                val = join(str,"");
                X(row,4) = (table.transit_walk_1(i) + table.transit_walk_2(i)) * table{i,val};

                if (table.transit_mode(i) == 'Subway (T)')
                    % Waiting times
                    str = ["Subway_wait_",num2str(wait)]; 
                    val = join(str,"");
                    X(row,5) = table{i,val};

                    % Riding times
                    str = ["Subway_IVTT_",num2str(transit_ride)]; 
                    val = join(str,"");
                    X(row,7) = table.Distance(i) * table{i,val};

                    % Price
                    X(row,10) = table.Transit_cost(i);

                elseif (table.transit_mode(i) == 'Bus')
                    % Waiting times
                    str = ["Bus_wait_",num2str(wait)]; 
                    val = join(str,"");
                    X(row,5) = table{i,val};

                    % Riding times
                    str = ["Bus_IVTT_",num2str(transit_ride)]; 
                    val = join(str,"");
                    X(row,7) = table.Distance(i) * table{i,val};

                    % Price
                    X(row,10) = table.Transit_cost(i);
                    
                else
                    % Waiting times
                    str = ["Rail_wait_",num2str(wait)]; 
                    val = join(str,"");
                    X(row,5) = table{i,val};

                    % Riding times
                    str = ["Rail_IVTT_",num2str(transit_ride)]; 
                    val = join(str,"");
                    X(row,7) = table.Distance(i) * table{i,val};

                    % Price
                    X(row,10) = table.Transit_cost(i);
                end
                
                X(row,12) = 1;
                
            elseif k == 3 % Exclusive rideshare alternative
                if(table{i,val1} == 'Exclusive rideshare')
                    X(row,3) = 1;
                end   
                
                % Walking times
                % Either 0 min (but this makes X singular) or drawn from a uniform random dist between 0-1 min
                X(row,4) = 0;

                % Waiting times
                str = ["Exc_wait_",num2str(wait)]; 
                val = join(str,"");
                X(row,5) = table{i,val};

                % Riding times
                str = ["Exc_IVTT_",num2str(exc_ride)]; 
                val = join(str,"");
                X(row,8) = table.Distance(i) * table{i,val};

                % Price
                str = ["Exc_cost_variability",num2str(exc_cost_var)]; 
                val = join(str,"");

                str2 = ["Exc_driver_wait",num2str(exc_driver_wait)]; 
                val2 = join(str2,"");
                X(row,10) = table{i,val} * (2.2 + table.Exc_cost(i) * table.Distance(i) + 0.42 * table{i,val2});
        
                X(row,13) = 1;
                
            elseif k == 4 % Pooled rideshare alternative
                if(table{i,val1} == 'Pooled rideshare')
                    X(row,3) = 1;
                end   
                
                % Walking times
                str = ["Pooled_walk_",num2str(walk)]; 
                val = join(str,"");
                X(row,4) = table{i,val};

                % Waiting times
                str = ["Pooled_wait_",num2str(wait)]; 
                val = join(str,"");
                X(row,5) = table{i,val};

                % Riding times
                str = ["Pooled_IVTT",num2str(pool_ride)]; 
                val = join(str,"");
                X(row,9) = table.Distance(i) * table{i,val};

                % Price
                str = ["Pool_cost_variability",num2str(pool_cost_var)]; 
                val = join(str,"");

                X(row,10) = table{i,val} * (2.2 + table.Pool_cost(i) * table.Distance(i));
                X(row,14) = 1;
            end
            row = row + 1;
        end 
        choice = choice + 1;
    end
end

writematrix(X,'Xmat.txt');
  
    