function pi = prob_weight(p,alpha,weight_type)
    if strcmp(weight_type,'original')
        pi = exp(-(-log(p)).^alpha);
    elseif strcmp(weight_type,'modified')
        pi = (p.^alpha)/(p.^alpha + (1-p).^alpha).^(1/alpha);
    end
end