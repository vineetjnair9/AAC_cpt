function [new_opt_soln,new_opt_cost,new_nu,rel_sens_soln,rel_sens_cost] = derivations_global(theta,nominal,scenario_vals,devs)
% theta: Parameter w.r.t which we compute sensitivity
% nominal: Nominal parameter values
% scenario: Travel times, prices, utility coefficients etc. for hypothetical choice situation
% active: Which constraints are active
% devs: Amount of deviations/perturbations in parameter value (as %)

    %% Set nominal and scenario values
    alpha_val = nominal(1);
    beta_val = nominal(2);
    lambda_val = nominal(3);
    p_val = nominal(4);    
    
    b_sm = scenario_vals.b_sm;
    R_val = scenario_vals.R;
    lb = scenario_vals.lb;
    ub = scenario_vals.ub;
    u0 = scenario_vals.u0;
    x1 = scenario_vals.x1; % u_sm = x_sm + b * gamma
    x2 = scenario_vals.x2;
    
    %% Calculate symbolic derivatives and expressions
    gamma = sym('gamma'); 
    alpha = sym('alpha'); 
    beta = sym('beta'); 
    lambda = sym('lambda'); 
    p = sym('p'); 
    nu1 = sym('nu1'); 
    nu2 = sym('nu2');

    u1 = x1 + b_sm*gamma;
    u2 = x2 + b_sm*gamma;
    
    % Type of reference used
    if R_val == 'Best case' % dynamic
        R = u2;
        U_sR = -exp(-(-log(p))^alpha) * lambda * (R - u1)^beta; 
        A_sR = -lambda*(R - u0)^beta;
    elseif R_val == 'Worst case' % dynamic
        R = u1;
        U_sR = exp(-(-log(p))^alpha) * (u2 - R)^beta;
        A_sR = (u0 - R)^beta;
    end
    
    p_sR = exp(U_sR)/(exp(U_sR) + exp(A_sR));
    f = gamma*p_sR;

    % Derivatives wrt price
    p_sR_dot = diff(p_sR,gamma); % w.r.t gamma
    p_sR_ddot = diff(p_sR_dot,gamma);
    
    theta_val = str2sym(theta);
    
    if theta == 'alpha'
        theta_nom = alpha_val;
        f = subs(f,[alpha, beta, lambda, p],[theta_val, beta_val, lambda_val, p_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[theta_val, beta_val, lambda_val, p_val]);
    elseif theta == 'beta'
        theta_nom = beta_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, theta_val, lambda_val, p_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, theta_val, lambda_val, p_val]);
    elseif theta == 'lambda'
        theta_nom = lambda_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, beta_val, theta_val, p_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, beta_val, theta_val, p_val]);
    else 
        theta_nom = p_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, beta_val, lambda_val, theta_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, beta_val, lambda_val, theta_val]);
    end
 
    f_dot = diff(f,gamma); % gamma*p_sR_dot + p_sR
    f_ddot = diff(f_dot,gamma);
           
    perturbations = length(devs);
    opt_check = zeros(1,perturbations);
    prob_actual = zeros(1,perturbations); % Using optimal price matching actual parameters    opt_cost_actual = zeros(1,length(devs));
    opt_cost = zeros(1,perturbations);
    opt_soln = zeros(1,perturbations);
    delta_soln = zeros(1,perturbations);
    delta_cost = zeros(1,perturbations);
    theta_vals = zeros(1,perturbations);
    delta_theta = zeros(1,perturbations);
    theta_noms = zeros(1,perturbations); % Keeps track of current nominal point - local domain around this point similar to basin of attraction
    rel_sens_soln = zeros(1,perturbations);
    rel_sens_cost = zeros(1,perturbations);
   
    [gamma_nom,f_nom_val,dGamma_dTheta_nom,dF_dTheta_nom,nu1,nu2,max_dev] = derivations_local_module(f,theta,theta_nom);
    
    for i = 1:1:perturbations
        if (devs(i) <= max_dev) % Can treat as local changes
            theta_noms(i) = theta_nom;
            delta_theta(i)= theta_nom * (devs(i)/100); % theta - theta_nom
            theta_vals(i) = theta_nom + delta_theta(i);
            opt_soln(i) = gamma_nom + (dGamma_dTheta_nom * delta_theta(i));
            delta_soln(i) = opt_soln(i) - gamma_nom;
            opt_check(i) = double(subs(L_gamma,[gamma,theta_val],[opt_soln(i),theta_vals(i)]));
            opt_cost(i) = subs(f,[gamma,theta_val],[opt_soln(i),theta_vals(i)]);
            prob_actual(i) = double(subs(p_sR,[gamma,theta_val],[opt_soln(i),theta_vals(i)]));
            delta_cost(i) = opt_cost(i) - f_nom_val;
            %% Relative sensitivities
            % Relative sensitivities will be equal/constant for all perturbations
            % if we use linear 1st order Taylor approximation
            rel_sens_soln(i) = (delta_soln(i)/delta_theta(i))*100;
            rel_sens_cost(i) = (delta_cost(i)/delta_theta(i)).*100;    

            %% Error in 1st order approximation of optimal solution 
            % For now just compute error by comparing with numerical solution?
            
        else % Active set changes
            
            [gamma_nom,f_nom_val,dGamma_dTheta_nom,dF_dTheta_nom,nu1,nu2,max_dev] = derivations_local_module(f,theta,theta_nom);

    %% Local sensitivity plots
    figure(1)
    subplot(2,2,1)
    % Plug optimal solution into actual objective function (always valid)
    plot(theta_vals,opt_cost_actual);
    xlabel('Parameter value');
    ylabel('Optimal value of objective function ($)')    
    
    subplot(2,2,2)
    plot(theta_vals,delta_cost_actual);
    xlabel('Parameter value');
    ylabel('Deviation from nominal optimum objective')
    
    subplot(2,2,3)
    plot(theta_vals,opt_soln);
    xlabel('Parameter value');
    ylabel('Optimal solution or tariff ($)')
    
    subplot(2,2,4)
    plot(theta_vals,delta_soln);
    xlabel('Parameter value');
    ylabel('Deviation from nominal optimum solution')
    
    figure(2)
    % Linearized version (valid only for small perturbations)   
    subplot(2,1,1)
    plot(theta_vals,opt_cost_linear);
    xlabel('Parameter value');
    ylabel('Optimal value of (linearized) objective function ($)')
    
    subplot(2,1,2)
    plot(theta_vals,delta_cost_linear);
    xlabel('Parameter value');
    ylabel('Deviation from nominal optimum (linearized) objective')
    
    figure(3)
    plot(theta_vals,prob_actual);
    xlabel('Parameter value');
    ylabel('Subjective probability of acceptance')
       
end
    

