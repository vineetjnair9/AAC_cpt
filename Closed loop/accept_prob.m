function p_sR = accept_prob(gamma,u0,x1,x2,p,b_sm,R_type,weight_type,alpha_plus,alpha_minus,beta_plus,beta_minus,lambda)
u1 = x1 + b_sm.*gamma;
u2 = x2 + b_sm.*gamma;
R = reference(u0,u1,u2,p,R_type);
U_sR = subjective(alpha_plus,alpha_minus,beta_plus,beta_minus,lambda,u1,u2,p,R,weight_type);
A_sR = value_func(beta_plus,beta_minus,lambda,u0,R);

p_sR = exp(U_sR)./(exp(U_sR) + exp(A_sR));
% p_sR = 1./(1 + exp(A_sR - U_sR));

end