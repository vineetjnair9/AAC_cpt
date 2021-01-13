function [xf, converged, x] = newtonNdGMRES_be(fhand,x0,p,u,dt,t)
% INPUTS:
%   fhand - function handle
%   x0 - initial guess
%   itpause - parameter for plotting
%   u - input function handle 
%
% Use Newton's method to solve the nonlinear function
% defined by function handle fhand with initial guess x0.
% Returns solution xf and all intermediate solutions
% x = [x1...xk...xf].
%
% itpause is parameter for plotting and defines the
% number of Newton steps that are plotted sequentially
% pauses between sub-steps.

tolf = 1e-5;         % function convergence tolerance 1e-10
tolx = 1e-5;          % step convergence tolerance 1e-8
maxIters = 10;       % max # of iterations
x00 = x0;             % initial guess

len = length(x0);

% Newton loop
for iter = 1:maxIters
    f = fhand(x0,u(t),t);               % evaluate function
    J = Jf_FiniteDifference('human_car_behaviour_v5',x0,p,u(t),t);
    J_be = eye(len) - dt*J;
%     [dx,~] = gmres(J_be,-f);
    dx = -J_be\f; 
    nf(iter) = norm(f,Inf);      % norm of f at step k+1
    ndx(iter) = norm(dx,Inf);    % norm of dx at step k+1
    x(:,iter) = x0 + dx;         % solution x at step k+1
    x0 = x(:,iter);              % set value for next guess
    if nf(iter) < tolf && ndx(iter) < tolx % check for convergence
        %fprintf('Converged in %d iterations\n',iter);
        converged = 1;
        break
    end
end
xf = x(:,iter);

if iter == maxIters % check for non-convergence
    fprintf('Non-Convergence after %d iterations!!!\n',iter);
    converged = 0;
end

x = [x00,x];
