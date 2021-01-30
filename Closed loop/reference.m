%% Determining reference point
function R = reference(u0,u1,u2,p,R_type)
     
    if  strcmp(R_type,'dynamic_SMoDS')
        % R = mean of 2 SMoDS outcomes
        R = p.*u1 + (1-p).*u2;
    elseif strcmp(R_type,'dynamic_u0') % Setting R as the certainty equivalent itself 
        % i.e. R = objective utility of current most frequent alternative
        R = u0;
%     elseif strcmp(R_type,'dynamic_mean_prob')
%         % R = average of current alternative and mean SMoDS outcome
%         % Accounting for probabilities
%         R = (u0 + p.*u1 + (1-p).*u2)./2;  
%     elseif strcmp(R_type,'dynamic_mean')
%         % R = simple mean of all 3 outcomes (no probabilities)
%         R = (u0 + u1 + u2)./3;
    elseif strcmp(R_type,'SMoDS_worst')
        R = u1;
    elseif strcmp(R_type,'SMoDS_best')
        R = u2;
    end
     
end