function [xf, converged, x] = newtonNd_trap_cont(fhand,x0,p,u,dt,t)
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

tolf = 1e-5;         % function convergence tolerance 
tolx = 1e-5;          % step convergence tolerance 
maxIters = 500;       % max # of iterations
x00 = x0;             % initial guess
dq = 0.001;

len = length(x0);
q = 0;

% Newton loop
while q <= 1
    for iter = 1:maxIters
        f = fhand(x0,u(t),t,q);               % evaluate function
        %J = jacobian_finite_difference('human_car_behaviour_v5',x0, p, 'constant_speed_input', t, epsilon);
        J = Jf_FiniteDifference('human_car_behaviour_v5',x0,p,u(t),t);
        J_trap = eye(len) - q*(dt/2)*J;
%         cond_no = cond(J_trap)
        %[dx,~] = gmres(J_trap,-f,);
        dx = -J_trap\f; 
        nf(iter) = norm(f,Inf);      % norm of f at step k+1
        ndx(iter) = norm(dx,Inf);    % norm of dx at step k+1
        x(:,iter) = x0 + dx;         % solution x at step k+1
        x0 = x(:,iter);              % set value for next guess
        q = q + dq;

        if nf(iter) < tolf && ndx(iter) < tolx % check for convergence
            %fprintf('Converged in %d iterations\n',iter);
            converged = 1;
            break
        end
    end
end
xf = x(:,iter);

if iter == maxIters % check for non-convergence
    fprintf('Non-Convergence after %d iterations!!!\n',iter);
    converged = 0;
end

x = [x00,x];
