function [x,converged,errf_k,errDeltax_k,relDeltax_k,iterations] = NewtonMethod(eval_f,x0,p,u,errf,errDeltax,relDeltax,MaxIter,visualize,FiniteDifference,eval_Jf)
% uses Newton Method to solve the nonlinear system f(x)=0
% x0 is the initial guess for Newton iteration
% p is a structure containing all parameters needed to evaluate f( )
% u is a vector with the values of all inputs at the time of interest
% eval_f = text string with name of function evaluating f for a given vector x 
% eval_Jf = text string with name of function evaluating Jacobian of f at vector x
% errF      = absolute equation error: how close do you want f to zero?
% errDeltax = absolute output error: how close do you want x?
% relDeltax = relative output error: how close do you want x in perentage?
% declares convergence if ALL three criteria are satisfied 
% MaxIter = maximum number of iterations allowed
% 
% [x,converged,errf_k,errDeltax_k,relDeltax_k,iterations] = NewtonMethod(eval_f,x0,p,u,errf,errDeltax,relDeltax,MaxIter,visualize,FiniteDifference,eval_Jf)

% copyright Luca Daniel, MIT 2018

k=1;
X(:,1) = x0;
f = feval(eval_f,x0,p,u);
Deltax = 0;
errf_k = norm(f,inf);
errDeltax_k = 0;
relDeltax_k = 0;
if visualize
   visualizeResults(1,X,1,'.b');
end
while k<MaxIter & (errf_k>errf | errDeltax_k>errDeltax | relDeltax_k>relDeltax),
   if FiniteDifference
      Jf = eval_Jf_FiniteDifference(eval_f,X(:,k),p,u);
   else 
      Jf = feval(eval_Jf,X(:,k),p,u);
   end
   Deltax = Jf\(-f);
   X(:,k+1)= X(:,k)+Deltax;
   k = k+1;
   if visualize
      visualizeResults([1:1:k],X,k,'.b');
   end
   f = feval(eval_f,X(:,k),p,u);
   errf_k = norm(f,inf);
   errDeltax_k = norm(Deltax,inf);
   relDeltax_k = norm(Deltax,inf)/max(abs(X(:,k)));
end

x = X(:,k);

if errf_k<=errf & errDeltax_k<=errDeltax & relDeltax_k<=relDeltax
   converged = 1;
else 
   converged = 0;
end

iterations = k-1;