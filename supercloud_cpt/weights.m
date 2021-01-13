%% Calculating subjective weights
function w = weights(p,u1,u2,R,alpha_plus,alpha_minus,weight_type,cdf)
    F = @(u) bernoulli_cdf(p,u,u1,u2);
    pi_gain = @(x) distort_gain(x,alpha_plus,weight_type);
    pi_loss = @(x) distort_loss(x,alpha_minus,weight_type);
    
    if (cdf)
        if (u1 < R)
            w(1) = pi_loss(F(u1));
        else
            w(1) = 1 - pi_gain(1-F(u1));
        end

        if (u2 < R) 
            w(2) = pi_loss(F(u2)) - pi_loss(F(u1));
        else
            w(2) = pi_gain(1-F(u1)) - pi_gain(1-F(u2));
        end
    else
        if (u1 < R)
            w(1) = pi_loss(p);
        else
            w(1) = pi_gain(p);
        end
    
        if (u2 < R)
            w(2) = pi_loss(1-p);
        else
            w(2) = pi_gain(1-p);
        end
    end
end
