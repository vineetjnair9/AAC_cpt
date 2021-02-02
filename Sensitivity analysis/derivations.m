%% CDC Math derivations

% Varying only one factor at a time

function [opt_soln, opt_multiplier, opt_cost] = derivations(theta,nominal,scenario_vals,active,dev)
% theta: Parameter w.r.t which we compute sensitivity
% nominal: Nominal parameter values
% scenario: Travel times, prices, utility coefficients etc. for hypothetical choice situation
% active: Which constraints are active
% dev: Amount of deviatio/perturbation in parameter value (as %)

    %% Set nominal and scenario values
    alpha_val = nominal(1);
    beta_val = nominal(2);
    lambda_val = nominal(3);
    p_val = nominal(4);    

    gamma_0_val = scenario_vals.gamma_0;
    a_sm_val = scenarios_vals.a_sm;
    a_0_val = scenarios_vals.a_0;
    b_sm_val = scenarios_vals.b_sm;
    b_0_val = scenarios_vals.b_0;
    c_sm_val = scenarios_vals.c_sm;
    c_0_val = scenarios_vals.c_0;
    t0_val = scenarios_vals.t0;
    t1_val = scenarios_vals.t1;
    t2_val = scenario_vals.t2;
    R_val = scenario_vals.R;
    lb_val = scenarios_vals.lb;
    ub_val = scenarios_vals.lb;
    
    %% Calculate symbolic derivatives and expressions
    syms gamma alpha beta lambda p gamma_0 a_sm a_0 b_sm b_0 c_sm c_0 t0 t1 t2 R lb ub nu1 nu2 

    u0 = a_0*t0 + b_0*gamma_0 + c_0;
    u1 = a_sm*t1 + b_sm*gamma + c_sm;
    u2 = a_sm*t2 + b_sm*gamma + c_sm;
    
    % Type of reference used
    if R_val == 'Best case'
        R = u2;
        U_sR = -exp(-(-log(p))^alpha) * lambda * (R - u1)^beta; 
        A_sR = -lambda*(R - u0)^beta;
    elseif R_val == 'Worst case'
        R = u1;
        % TODO: Need to finish case for all possible static/dynamic reference types!
    end
    
    p_sR = exp(U_sR)/(exp(U_sR) + exp(A_sR));
    f = gamma*p_sR;

    % Derivatives wrt price
    p_sR_dot = diff(p_sR,gamma); % w.r.t gamma
    p_sR_ddot = diff(p_sR_dot,gamma);
    
    theta_val = str2sym(theta);
    
    if theta == 'alpha'
        theta_nom = alpha_val;
        f = subs(f,[alpha, beta, lambda, p, gamma_0, a_sm, a_0, b_sm, b_0, c_sm, c_0, t0, t1, t2],...
            [theta_val, beta_val, lambda_val, p_val, gamma_0_val, a_sm_val, a_0_val, b_sm_val,...
            b_0_val, c_sm_val, c_0_val, t0_val, t1_val, t2_val]);
        f_nom = subs(f,alpha,theta_nom);
    elseif theta == 'beta'
        theta_nom = beta_val;
        f = subs(f,[alpha, beta, lambda, p, gamma_0, a_sm, a_0, b_sm, b_0, c_sm, c_0, t0, t1, t2],...
            [alpha_val, theta_val, lambda_val, p_val, gamma_0_val, a_sm_val, a_0_val, b_sm_val,...
            b_0_val, c_sm_val, c_0_val, t0_val, t1_val, t2_val]);
        f_nom = subs(f,beta,theta_nom);
    elseif theta == 'lambda'
        theta_nom = lambda_val;
        f = subs(f,[alpha, beta, lambda, p, gamma_0, a_sm, a_0, b_sm, b_0, c_sm, c_0, t0, t1, t2],...
            [alpha_val, beta_val, theta_val, p_val, gamma_0_val, a_sm_val, a_0_val, b_sm_val,...
            b_0_val, c_sm_val, c_0_val, t0_val, t1_val, t2_val]);
        f_nom = subs(f,lambda,theta_nom);
    else 
        theta_nom = p_val;
        f = subs(f,[alpha, beta, lambda, p, gamma_0, a_sm, a_0, b_sm, b_0, c_sm, c_0, t0, t1, t2],...
            [alpha_val, beta_val, lambda_val, theta_val, gamma_0_val, a_sm_val, a_0_val, b_sm_val,...
            b_0_val, c_sm_val, c_0_val, t0_val, t1_val, t2_val]);
        f_nom = subs(f,p,theta_nom);
    end
 
    f_dot = diff(f,gamma); % gamma*p_sR_dot + p_sR
    f_ddot = diff(f_dot,gamma);

    % TODO: Need to check if this Hessian of the Lagrangian is indeed PD for SSC to hold!
    
    %% Calculate nominal optimum    
    % Can use CVX if -f is convex i.e. f is concave
    
    cvx_begin
        variable gamma_nom(1) % Nominal optimal solution
        f_nom = subs(f_nom,gamma,gamma_nom);
        minimize -f_0    
        subject to
            gamma_nom <= ub_val;
            gamma_nom >= lb_val;
    cvx_end
    
    f_nom = cvx_optval; % Nominal value of objective function
    
    %% Need to use other solver if problem not convex
    % Using MATLAB's inbuilt solver
    [gamma_nom, f_nom, exitflag] = fmincon(f_nom,(lb_val + ub_val)/2,[1; -1],[ub_val, -lb_val]);
    
    %% Local sensitivity analysis - assuming active set doesn't change
    if active == 1
        % Active constraint = upper bound (nu1 > 0, nu2 = 0)
        L = f + nu1*(gamma-ub); % Lagrangian
        L_gamma = diff(L,gamma);
        L_gamma2 = diff(L_gamma,gamma);
        
        KKT_matrix = [L_gamma2 1;1 0];
        L_gamma_theta = diff(L_gamma,theta_val);
        RHS = [L_gamma_theta; 0];

        % Sensitivity differentials (of optimal solution & Lagrange multipliers)
        result = -inv(KKT_matrix) * RHS;
        dGamma_dTheta = result(1);
        dLambda_dTheta = result(2);
        
    elseif active == -1 
        % Active constraint = lower bound (nu1 = 0, nu2 > 0)
        L = f + nu2*(lb-gamma); % Lagrangian
        L_gamma = diff(L,gamma);
        L_gamma2 = diff(L_gamma,gamma);
        
        KKT_matrix = [L_gamma2 -1;-1 0];
        L_gamma_theta = diff(L_gamma,theta_val);
        RHS = [L_gamma_theta; 0];

        % Sensitivity differentials (of optimal solution & Lagrange multipliers)
        result = -inv(KKT_matrix) * RHS;
        dGamma_dTheta = result(1);
        dLambda_dTheta = result(2);
        
    else 
        % neither constraint is active
        L = f;
        L_gamma = diff(L,gamma);
        L_gamma2 = diff(L_gamma,gamma);
        L_gamma_theta = diff(L_gamma,theta_val);
        dGamma_dTheta = -inv(L_gamma2) * L_gamma_theta;
    end
   
    % Evaluate differentials at the optimum point
    gGamma_dTheta_nom = subs(dGamma_dTheta,[gamma,theta_val],[gamma_nom,theta_nom]);

    delta_theta = theta_nom * (dev/100);
    opt_soln = gamma_nom + (gGamma_dTheta_nom * delta_theta);
    
    
    
end


    

