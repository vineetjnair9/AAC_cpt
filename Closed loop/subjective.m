%% Calculating subjective utilities
function U_sR = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p1,R,weight_type) 

    V = @(u) value_func(beta_plus,beta_minus,lambda,u,R);
    u = [u1;u2];
    p_vals = [p1,1-p1];

    [u_low,ind_low] = min(u,[],1);
    [u_high,ind_high] = max(u,[],1);

    p = p_vals(ind_low);

    [w1,w2] = weights(p,u_low,u_high,R,alpha_plus,alpha_minus,weight_type);
    U_sR = w1 .* V(u_low) + w2 .* V(u_high);
end
   