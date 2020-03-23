%% CDC Math derivations

% Varying only one factor at a time
% V2 - Using lumped up versions for objective utilities 
% Also modified symbolic vs numeric variables

function [gamma_nom,f_nom_val,dGamma_dTheta_nom,dF_dTheta_nom,nu1,nu2,max_dev] = derivations_local_module(f,theta,theta_nom)
% theta: Parameter w.r.t which we compute sensitivity
% nominal: Nominal parameter values
% scenario: Travel times, prices, utility coefficients etc. for hypothetical choice situation
% active: Which constraints are active
% devs: Amount of deviations/perturbations in parameter value (as %)

    if theta == 'alpha'
        f_nom = subs(f,alpha,theta_nom);
    elseif theta == 'beta'
        f_nom = subs(f,beta,theta_nom);
    elseif theta == 'lambda'
        f_nom = subs(f,lambda,theta_nom);
    else 
        f_nom = subs(f,p,theta_nom);
    end
 
    f_dot = diff(f,gamma); % gamma*p_sR_dot + p_sR
    f_ddot = diff(f_dot,gamma);
   
    %% Using MATLAB's inbuilt solver
    gamma_ = sym('gamma_');
    f_nom = subs(f_nom,gamma,gamma_);
    f_nom = matlabFunction(f_nom);
    minusf_nom = @(gamma_) -f_nom(gamma_);
    options = optimoptions(@fmincon, 'Algorithm','interior-point');
    
    init = (lb + ub)/2;
    
    % Global minimum
    gs = GlobalSearch;
    
    problem = createOptimProblem('fmincon','objective',...
    minusf_nom,'x0',init,'lb',lb,'ub',ub,'options',options);
    
    % Plot out objective function f_nom to visualize local & global minima
%     gamma_vals = linspace(lb,ub,100);
%     plot(gamma_vals,f_nom(gamma_vals))
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

end
    

