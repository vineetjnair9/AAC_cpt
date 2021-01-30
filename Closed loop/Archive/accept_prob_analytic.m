function [p_sR] = accept_prob_analytic(gamma,alpha,beta,lambda,x1,x2,p,b_sm,u0)
% p_sR = 1./(1 + exp(lambda.*(exp(-(-log(p)).^alpha).*(x2 - x1).^beta - (x2 + b_sm.*gamma - u0).^beta)));
p_sR = 1./(1 + exp(lambda.*(exp(-(-log(p)).^alpha).*(x2 - x1).^beta - abs(x2 + b_sm.*gamma - u0).^beta)));
end