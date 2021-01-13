%% CDC Math derivations
% Varying only one factor at a time 
% Simply resolves whole optimization problem at each new operating point
% for new solution & new cost/objective function value

function [gamma_opt,f_opt,rel_sens_soln,rel_sens_cost,mismatch_loss,theta_vals,theta_nom] = derivations_numerical(theta,nominal,scenario_vals,devs,plots)
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
    if strcmp(R_val,'Best case')
        R = u2;
        U_sR = -exp(-(-log(p))^alpha)*lambda*(R - u1)^beta; 
        A_sR = -lambda*(R - u0)^beta;
        
    elseif strcmp(R_val,'Worst case')
        R = u1;
        U_sR = exp(-(-log(p))^alpha)*(u2 - R)^beta;
        A_sR = (u0 - R)^beta;
               
    elseif strcmp(R_val,'Alternative')
        R = u0;
        A_sR = 0;
        U_sR = -exp(-(-log(p))^alpha)*lambda*(R - u1)^beta + exp(-(-log(1-p))^alpha) * (u2 - R)^beta;       
    
    elseif strcmp(R_val,'Mean')
        R = p_val*u1 + (1-p)*u2;
        U_sR = -exp(-(-log(p))^alpha)*lambda*(R - u1)^beta + exp(-(-log(1-p))^alpha) * (u2 - R)^beta;
        if (u0 < R)      
            A_sR = -lambda*(R - u0)^beta;
        else
            A_sR = (u0 - R)^beta;
        end
    end   
    
    p_sR = exp(U_sR)/(exp(U_sR) + exp(A_sR));
    f = gamma*p_sR;

    % Derivatives wrt price
    p_sR_dot = diff(p_sR,gamma); % w.r.t gamma
    p_sR_ddot = diff(p_sR_dot,gamma);
    
    theta_val = str2sym(theta);
    
    if strcmp(theta,'alpha')
        theta_nom = alpha_val;
        f = subs(f,[alpha, beta, lambda, p],[theta_val, beta_val, lambda_val, p_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[theta_val, beta_val, lambda_val, p_val]);
    elseif strcmp(theta,'beta')
        theta_nom = beta_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, theta_val, lambda_val, p_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, theta_val, lambda_val, p_val]);
    elseif strcmp(theta,'lambda')
        theta_nom = lambda_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, beta_val, theta_val, p_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, beta_val, theta_val, p_val]);
    else 
        theta_nom = p_val;
        f = subs(f,[alpha, beta, lambda, p],[alpha_val, beta_val, lambda_val, theta_val]);
        p_sR = subs(p_sR,[alpha, beta, lambda, p],[alpha_val, beta_val, lambda_val, theta_val]);
    end

    %% Calculate nominal optimum     
    gamma_ = sym('gamma_');
    
    options = optimoptions(@fmincon, 'Algorithm','interior-point');
    init = (lb + ub)/2;
    % Global minimum
    gs = GlobalSearch;
        
    f_nom = subs(f,[gamma,theta_val],[gamma_,theta_nom]);
    f_nom = matlabFunction(f_nom); % Expected revenue
    minusf_nom = @(gamma_) -f_nom(gamma_);
    problem = createOptimProblem('fmincon','objective',...
            minusf_nom,'x0',init,'lb',lb,'ub',ub,'options',options);
    [gamma_nom, f_nom_val] = run(gs,problem);

    f_nom_val = -f_nom_val;
    delta_theta = theta_nom .* (devs./100); % theta - theta_nom
    theta_vals = theta_nom + delta_theta;
    perturbations = length(devs);
    gamma_opt = zeros(1,perturbations);
    f_opt = zeros(1,perturbations);
    prob_actual = zeros(1,perturbations);
    mismatch_loss = zeros(1,perturbations);
    
    for i = 1:perturbations
        f_obj = subs(f,[gamma,theta_val],[gamma_,theta_vals(i)]);
        f_obj = matlabFunction(f_obj); % Expected revenue
        minusf_obj = @(gamma_) -f_obj(gamma_);
        problem = createOptimProblem('fmincon','objective',...
            minusf_obj,'x0',init,'lb',lb,'ub',ub,'options',options);
        [gamma_opt(i), f_opt(i)] = run(gs,problem);
        prob_actual(i) = double(subs(p_sR,[gamma,theta_val],[gamma_opt(i),theta_vals(i)]));
        mismatch_loss(i) = f_nom_val - double(subs(f,[gamma,theta_val],[gamma_opt(i),theta_nom]));
    end
    
    f_opt = -f_opt;
    
    delta_soln = gamma_opt - gamma_nom;
    delta_cost = f_opt - f_nom_val;
    
    %% Relative sensitivities
    % Relative sensitivities will be equal/constant for all perturbations
    % if we use linear 1st order Taylor approximation
    rel_sens_soln = (delta_soln./delta_theta).*100;
    rel_sens_cost = (delta_cost./delta_theta).*100;    
    
    %%  Plots
    if plots == 1
        figure(1)
        subplot(2,1,1)
        % Plug optimal solution into actual objective function (always valid)
        plot(theta_vals,f_opt);
        xlabel('Parameter value');
        ylabel('Optimal value of objective function ($)')    

        subplot(2,1,2)
        plot(theta_vals,gamma_opt);
        xlabel('Parameter value');
        ylabel('Optimal solution or tariff ($)')
        figure(2)

        figure(2)
        plot(theta_vals,prob_actual);
        xlabel('Parameter value');
        ylabel('Subjective probability of acceptance')
    end      
end