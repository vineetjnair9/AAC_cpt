%% PARAMETER ESTIMATION
% Vineet Jagadeesan Nair

% Maximum simulated likelihood estimation by Kenneth Train
% Assumptions:
% (1) All coefficients as random
% (2) Using Bernoulli CDF

%% I - Estimation of mode-choice model
% Multinomial logit for utility functions
% Generation of data matrix X for using Kenneth Train's code

clear; clc;
T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_6-10_n955');

load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Code\AAC_cpt\Survey analysis\factor_levels.mat')

%%
table_orig = T.Fulllaunch610n955;

% Vary if table columns change
table = convertvars(table_orig,[7:8,10:12,48:139],'double');

num_responses = size(table,1); % Total no. of survey responses
scenarios = 11*num_responses; % 11 choice scenarios per respondent
alternatives = scenarios*3; % 3 alternatives per choice situat ion

X = zeros(alternatives,12); % Observations for each scenario on 4 predictor variables            
option = 1;
choice = 1;

% Loop through each respondent
for i = 1:num_responses
    
    % Loop through each choice scenario
    for j = 1:11
        
        if table.transit_mode(i) == 'Bus'    
            str1 = ["bus",num2str(j)]; 
            val1 = join(str1,"");
        elseif table.transit_mode(i) == 'Subway (T)'
            str1 = ["subway",num2str(j)]; 
            val1 = join(str1,"");
        elseif table.transit_mode(i) == 'Commuter rail'
            str1 = ["Rail",num2str(j)]; 
            val1 = join(str1,"");
        end    
        
        % Find correct factor levels
        exc_cost_var = factor_levels(j,1);
        pool_cost_var = factor_levels(j,2);
        wait = factor_levels(j,3);
        walk = factor_levels(j,4);
        transit_ride = factor_levels(j,5);
        exc_ride = factor_levels(j,6);
        pool_ride = factor_levels(j,7);
        exc_driver_wait = factor_levels(j,8);       
              
        % Loop through each alternative
        for k = 1:3
            X(option,1) = i; % Person facing this alternative
            X(option,2) = choice; % Choice situations numbered sequentially            
        
            if k == 1 % Public transit alternative
                if(table{i,val1} == 'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}')
                    X(option,3) = 1;
                end    
                
                % Walking times
                str = ["Transit_walk_variability_",num2str(walk)]; 
                val = join(str,"");
                X(option,4) = (table.transit_walk_1(i) + table.transit_walk_2(i)) * table{i,val};

                if (table.transit_mode(i) == 'Subway (T)')
                    % Waiting times
                    str = ["Subway_wait_",num2str(wait)]; 
                    val = join(str,"");
                    X(option,5) = table{i,val};

                    % Riding times
                    str = ["Subway_IVTT_",num2str(transit_ride)]; 
                    val = join(str,"");
                    X(option,6) = table.Distance(i) * table{i,val};

                    % Price
                    X(option,9) = table.Transit_cost(i);

                elseif (table.transit_mode(i) == 'Bus')
                    % Waiting times
                    str = ["Bus_wait_",num2str(wait)]; 
                    val = join(str,"");
                    X(option,5) = table{i,val};

                    % Riding times
                    str = ["Bus_IVTT_",num2str(transit_ride)]; 
                    val = join(str,"");
                    X(option,6) = table.Distance(i) * table{i,val};

                    % Price
                    X(option,9) = table.Transit_cost(i);
                    
                else
                    % Waiting times
                    str = ["Rail_wait_",num2str(wait)]; 
                    val = join(str,"");
                    X(option,5) = table{i,val};

                    % Riding times
                    str = ["Rail_IVTT_",num2str(transit_ride)]; 
                    val = join(str,"");
                    X(option,6) = table.Distance(i) * table{i,val};

                    % Price
                    X(option,9) = table.Transit_cost(i);
                end
                
                X(option,10) = 1;
                
            elseif k == 2 % Exclusive rideshare alternative
                if(table{i,val1} == 'Exclusive rideshare')
                    X(option,3) = 1;
                end   
                
                % Walking times
                X(option,4) = 0;

                % Waiting times
                str = ["Exc_wait_",num2str(wait)]; 
                val = join(str,"");
                X(option,5) = table{i,val};

                % Riding times
                str = ["Exc_IVTT_",num2str(exc_ride)]; 
                val = join(str,"");
                X(option,7) = table.Distance(i) * table{i,val};

                % Price
                str = ["Exc_cost_variability",num2str(exc_cost_var)]; 
                val = join(str,"");

                str2 = ["Exc_driver_wait",num2str(exc_driver_wait)]; 
                val2 = join(str2,"");
                X(option,9) = table{i,val} * (2.2 + table.Exc_cost(i) * table.Distance(i) + 0.42 * table{i,val2});
        
                X(option,11) = 1;
                
            elseif k == 3 % Pooled rideshare alternative
                if(table{i,val1} == 'Pooled rideshare')
                    X(option,3) = 1;
                end   
                
                % Walking times
                str = ["Pooled_walk_",num2str(walk)]; 
                val = join(str,"");
                X(option,4) = table{i,val};

                % Waiting times
                str = ["Pooled_wait_",num2str(wait)]; 
                val = join(str,"");
                X(option,5) = table{i,val};

                % Riding times
                str = ["Pooled_IVTT",num2str(pool_ride)]; 
                val = join(str,"");
                X(option,8) = table.Distance(i) * table{i,val};

                % Price
                str = ["Pool_cost_variability",num2str(pool_cost_var)]; 
                val = join(str,"");

                X(option,9) = table{i,val} * (2.2 + table.Pool_cost(i) * table.Distance(i));
                X(option,12) = 1;
            end
            option = option + 1;
        end 
        choice = choice + 1;
    end
end

writematrix(X,'Xmat.txt');
  
    