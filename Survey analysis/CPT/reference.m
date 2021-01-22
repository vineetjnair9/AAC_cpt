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
        if strcmp(table.Reference(i),'Public transit: ${q://QID56/ChoiceGroup/SelectedChoices}')
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