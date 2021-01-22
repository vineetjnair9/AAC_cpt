%% Calculate error for nonlinear equation
function error = error_func(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,CE_actual,u1,u2,p,R,weight_type,cdf)
    CE_pred = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p,R,weight_type,cdf);
    
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
