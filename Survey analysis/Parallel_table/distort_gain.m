%% Probability weighting function
function pi = distort_gain(p,alpha_plus,weight_type)
    if strcmp(weight_type,'original')
        pi = exp(-(-log(p))^alpha_plus);
    elseif strcmp(weight_type,'modified')
        pi = (p^alpha_plus)/(p^alpha_plus + (1-p)^alpha_plus)^(1/alpha_plus);
    end
end

