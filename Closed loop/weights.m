%% Calculating subjective probability weights
function [w1,w2] = weights(p,u1,u2,R,alpha_plus,alpha_minus,weight_type)
    F = @(u) bernoulli_cdf(p,u,u1,u2);
    pi_gain = @(x) distort_gain(x,alpha_plus,weight_type);
    pi_loss = @(x) distort_loss(x,alpha_minus,weight_type);
    
    if (u1 < R)
        w1 = pi_loss(F(u1));
    else
        w1 = 1 - pi_gain(1-F(u1));
    end
    
    if (u2 < R) 
        w2 = pi_loss(F(u2)) - pi_loss(F(u1));
    else
        w2 = pi_gain(1-F(u1)) - pi_gain(1-F(u2));
    end
end