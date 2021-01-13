%% Calculating subjective utilities
function U_sR = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p1,R,weight_type,cdf) 

    V = @(u) value_func(beta_plus,beta_minus,lambda,u,R);
    u = [u1,u2];
    p_vals = [p1,1-p1];

    [u_low,ind_low] = min(u);
    [u_high,~] = max(u);

    p = p_vals(ind_low);

    w = weights(p,u_low,u_high,R,alpha_plus,alpha_minus,weight_type,cdf);
    U_sR = w(1) * V(u_low) + w(2) * V(u_high);
end
