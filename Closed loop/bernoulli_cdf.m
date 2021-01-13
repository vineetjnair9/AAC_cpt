%% Bernoulli CDF (binomial distribution with only 1 trial & 2 outcomes)
function F = bernoulli_cdf(p,u,u1,u2)

    a = u < u1;
    b = u >= u1 & u < u2;
    c = u > u2;
    F = zeros(size(u));
    F(a) = 0;
    F(b) = p(b);
    F(c) = 1;
    
%     if (u < u1)
%         F = 0;
%     elseif (u >= u1 & u < u2)
%         F = p;
%     else
%         F = 1;
%     end

end