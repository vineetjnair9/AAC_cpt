%% CDC Math derivations

% Varying only one factor at a time
% V2 - Using lumped up versions for objective utilities 
% Also modified symbolic vs numeric variables

function [dGamma_dTheta_nom,dF_dTheta_nom,opt_soln,opt_cost_actual,error,max_dev,rel_sens_soln,rel_sens_cost_actual] = derivations_local(theta,nominal,scenario_vals,devs,plots)
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
        
        % TODO: Need to finish cases for all possible static/dynamic reference types!
        % R_val = 'Mean/expected value' (dynamic)
        % R_val = 'Alternative' (static)
    end
    
%     if (u0 < R)      
%         A_sR = -lambda*(R - u0)^beta;
%     else
%         A_sR = (u0 - R)^beta;
%     end
    
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
        f_nom = subs(f,alpha,theta_nom);
    elseif theta == 'beta'
        theta_nom = beta_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, theta_val, lambda_val, p_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, theta_val, lambda_val, p_val]);
        f_nom = subs(f,beta,theta_nom);
    elseif theta == 'lambda'
        theta_nom = lambda_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, beta_val, theta_val, p_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, beta_val, theta_val, p_val]);
        f_nom = subs(f,lambda,theta_nom);
    else 
        theta_nom = p_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, beta_val, lambda_val, theta_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, beta_val, lambda_val, theta_val]);
        f_nom = subs(f,p,theta_nom);
    end
 
    f_dot = diff(f,gamma); % gamma*p_sR_dot + p_sR
    f_ddot = diff(f_dot,gamma);
    
    %% Calculate nominal optimum    
%     % Ensure that I find the global minimum (rather than one of many possible local minima)
%     
%     % Can use CVX if -f is convex i.e. f is concave
%     A = [1; -1];
%     b = [ub; -lb];
%     
%     gamma_nom = sym('gamma_nom');
%     f_nom = subs(f_nom,gamma,gamma_nom);
%     f_nom = matlabFunction(f_nom);
%     
%     cvx_begin
%         variable gamma_nom(1) % Nominal optimal solution
%         dual variable nu; % Lagrange multipliers        
%         maximize f_nom(gamma_nom)    
%         subject to
%             nu : A * gamma_nom <= b;
%     cvx_end
%     
%     nu1 = nu(1); nu2 = nu(2);
%     f_nom_val = cvx_optval; % Nominal value of objective function   
   
    %% Using MATLAB's inbuilt solver
    gamma_ = sym('gamma_');
    f_nom = subs(f_nom,gamma,gamma_);
    f_nom = matlabFunction(f_nom);
    minusf_nom = @(gamma_) -f_nom(gamma_);
    options = optimoptions(@fmincon, 'Algorithm','interior-point');
    
    init = (lb + ub)/2;
    
%     A = [1;-1]; 
%     b = [ub,-lb];
%     [gamma_nom, f_nom_val, exitflag,output,nu] = fmincon(f_nom,init,A,b);
%     nu1 = nu.ineqlin(1); nu2 = nu.ineqlin(2);
    
    % Local minimum
%     [gamma_nom, f_nom_val, exitflag,output,nu] = fmincon(-f_nom,init,[],[],[],[],lb,ub);
%     nu1 = round(nu.upper,7); 
%     nu2 = round(nu.lower,7);    
    
    % Global minimum
    gs = GlobalSearch;
    
    problem = createOptimProblem('fmincon','objective',...
    minusf_nom,'x0',init,'lb',lb,'ub',ub,'options',options);
    
    % Plot out objective function f_nom to visualize local & global minima
%     gamma_vals = linspace(lb,ub,100);
%     plot(gamma_vals,f_nom(gamma_vals));
%     xlabel('Price ($)'); ylabel('Expected revenue');

    [gamma_nom, f_nom_val, exitflag, output, manymins] = run(gs,problem);
    
    % Post-optimal calculation of Lagrange multiplers
    if (gamma_nom == ub) % upper bound active
        nu1 = subs(f_dot,[gamma,theta_val],[gamma_nom,theta_nom]);
        nu2 = 0;
    elseif (gamma_nom == lb)
        nu1 = 0;
        nu2 = -subs(f_dot,[gamma,theta_val],[gamma_nom,theta_nom]);
    else
        nu1 = 0;
        nu2 = 0;
    end
    
    %% Local sensitivity analysis - assuming active set doesn't change
     if nu1 > 0
        % Active constraint = upper bound (nu1 > 0, nu2 = 0)
        L = -f + nu1*(gamma-ub); % Lagrangian
        L_gamma = diff(L,gamma);
        L_gamma2 = diff(L_gamma,gamma);
        L_theta = diff(L,theta_val);
        
        L_gamma2_nom = double(subs(L_gamma2,[gamma,theta_val],[gamma_nom,theta_nom]));
        KKT_matrix = [L_gamma2_nom 1;1 0];
        
        L_gamma_theta = diff(L_gamma,theta_val);
        L_gamma_theta_nom = double(subs(L_gamma_theta,[gamma,theta_val],[gamma_nom,theta_nom]));
        RHS = [L_gamma_theta_nom; 0];

        % Sensitivity differentials (of optimal solution & Lagrange multipliers)
        result = -inv(KKT_matrix) * RHS;
        dGamma_dTheta_nom = result(1);
        dLambda_dTheta_nom = result(2);
        
        % Prediction of local sensitivity domain
        % max_dev = theta - theta_nom
        max_dev = -nu1/dLambda_dTheta_nom; % Linearization only valid for small perturbations
        max_dev = (max_dev/theta_nom)*100;
        
    elseif nu2 > 0
        % Active constraint = lower bound (nu1 = 0, nu2 > 0)
        L = -f + nu2*(lb-gamma); % Lagrangian
        L_gamma = diff(L,gamma);
        L_gamma2 = diff(L_gamma,gamma);
        L_theta = diff(L,theta_val);
        
        L_gamma2_nom = double(subs(L_gamma2,[gamma,theta_val],[gamma_nom,theta_nom]));
        KKT_matrix = [L_gamma2_nom -1;-1 0];
        
        L_gamma_theta = diff(L_gamma,theta_val);
        L_gamma_theta_nom = double(subs(L_gamma_theta,[gamma,theta_val],[gamma_nom,theta_nom]));
        RHS = [L_gamma_theta_nom; 0];

        % Sensitivity differentials (of optimal solution & Lagrange multipliers)
        result = -inv(KKT_matrix) * RHS;
        dGamma_dTheta_nom = result(1);
        dLambda_dTheta_nom = result(2);
        
        % Prediction of local sensitivity domain
        max_dev = -nu2/dLambda_dTheta_nom;
        max_dev = (max_dev/theta_nom)*100;
        
    else 
        % neither constraint is active (nu1 = nu2 = 0)
        L = -f;
        L_gamma = diff(L,gamma);
        L_gamma2 = diff(L_gamma,gamma);
        L_theta = diff(L,theta_val);
        L_gamma_theta = diff(L_gamma,theta_val);
        dGamma_dTheta = -inv(L_gamma2) * L_gamma_theta;
        
        % Prediction of local sensitivity domain
        % Both constraints are inactive - one of them enters active set
        
        dGamma_dTheta_nom = double(subs(dGamma_dTheta,[gamma,theta_val],[gamma_nom,theta_nom]));
        
        % If upper bound becomes active
        max_devs = zeros(1,2);
        max_devs(1) = -(gamma_nom - ub)/(dGamma_dTheta_nom);
        max_devs(1) = (max_devs(1)/theta_nom)*100;
        
        % If lower bound becomes active
        max_devs(2) = (lb - gamma_nom)/(dGamma_dTheta_nom);
        max_devs(2) = (max_devs(2)/theta_nom)*100;
        
        % max_dev = [min(max_devs),max(max_devs)];
        
        % Most conservative estimate of doman
        [~,index] = min(abs(max_devs));
        max_dev = max_devs(index);
      
    end

    %% Evaluate differentials at the optimum point
    dF_dTheta_nom = double(subs(L_theta,[gamma,theta_val],[gamma_nom,theta_nom]));
    L_gamma2_nom = double(subs(L_gamma2,[gamma,theta_val],[gamma_nom,theta_nom]));
    L_gamma_theta_nom = double(subs(L_gamma_theta,[gamma,theta_val],[gamma_nom,theta_nom]));
    L_theta2 = diff(L_theta,theta_val);
    L_theta2_nom = double(subs(L_theta2,[gamma,theta_val],[gamma_nom,theta_nom]));
    d2F_dTheta2_nom = L_gamma2_nom * (dGamma_dTheta_nom^2) + 2*L_gamma_theta_nom*dGamma_dTheta_nom + L_theta2_nom;
    
    %% 1st order sensitivity of optimal solution & optimal cost/value function
    delta_theta = theta_nom .* (devs./100); % theta - theta_nom
    theta_vals = theta_nom + delta_theta;
    
    opt_soln = gamma_nom + (dGamma_dTheta_nom .* delta_theta);
    delta_soln = opt_soln - gamma_nom;
    
    % Using formula from Buskens et al. 2000
    opt_cost_linear = f_nom_val + (dF_dTheta_nom .* delta_theta);
    delta_cost_linear = opt_cost_linear - f_nom_val;
    
    %% Check if 1st order optimality conditions still hold
    perturbations = length(devs);
    opt_check = zeros(1,perturbations);
    prob_actual = zeros(1,perturbations); % Using optimal price matching actual parameters    opt_cost_actual = zeros(1,length(devs));
    opt_cost_actual = zeros(1,perturbations);
    for i = 1:perturbations
        opt_price = opt_soln(i);
        param_value = theta_vals(i);
        opt_check(i) = double(subs(L_gamma,[gamma,theta_val],[opt_price,param_value]));
        opt_cost_actual(i) = subs(f,[gamma,theta_val],[opt_price,param_value]);
        prob_actual(i) = double(subs(p_sR,[gamma,theta_val],[opt_price,param_value]));
    end
    
    delta_cost_actual = opt_cost_actual - f_nom_val;
   
    %% Relative sensitivities
    % Relative sensitivities will be equal/constant for all perturbations
    % if we use linear 1st order Taylor approximation
    rel_sens_soln = (delta_soln./delta_theta).*100;
    rel_sens_cost_linear = (delta_cost_linear./delta_theta).*100;
    rel_sens_cost_actual = (delta_cost_actual./delta_theta).*100;    
    
    %% Estimate error in the 1st order approximation using 2nd order expansion 
    opt_cost2 = f_nom_val + (dGamma_dTheta_nom .* delta_theta) + 0.5 * d2F_dTheta2_nom .* (delta_theta).^2;
    error = ((opt_cost_linear - opt_cost2)./opt_cost2).*100; % Relative error (%) using 2nd order estimate as true value
   
    %% Local sensitivity plots
    if plots == 1
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
end
    

