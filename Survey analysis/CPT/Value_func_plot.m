%% Plotting value function
R = 0;
u = linspace(-100,100,10000);
p = linspace(0,1,10000);
V = zeros(1,length(u));
pi = zeros(1,length(p));

for i = 1:10000
    V(i) = value_fun1(beta,gamma,lambda,R,u(i));
    pi1(i) = exp(-(-log(p(i)))^alpha);
    pi2(i) 
end

figure(1)
plot(u,V)
ylabel('Subjective value V');
xlabel('Objective utility');
title('Value function')

figure(2)
plot(p,pi1)
ylabel('Subjective probability distortion');
xlabel('Objective probabilty');
title('Probability weighting function')
legend('riginal

%%
function V = value_fun1(beta,gamma,lambda,R,u)

if (u >= R)
    V = (u-R)^beta;
else
    V = -lambda*(R-u)^gamma;
end
end
