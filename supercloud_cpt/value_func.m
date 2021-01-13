%% Value function
function V = value_func(beta_plus,beta_minus,lambda,u,R)
    if (u >= R)
        V = (u-R).^beta_plus;
    else
        V = -lambda*(R-u).^beta_minus;
    end
end
