function pi = distort_loss(p,alpha_minus,weight_type)
    if strcmp(weight_type,'original')
        pi = exp(-(-log(p)).^alpha_minus);
    elseif strcmp(weight_type,'modified')
        pi = (p.^alpha_minus)/(p.^alpha_minus + (1-p).^alpha_minus).^(1/alpha_minus);
    end
end