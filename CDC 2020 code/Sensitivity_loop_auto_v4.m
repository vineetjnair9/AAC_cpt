%% Sensitivity analysis loop code v4

% Obtain both analytical and numerical results for optimal solution, objective and mismatch loss

% But randomize all scenario values like earlier
% Compare sensitivities w.r.t each of the 4 parameters
% Compare max allowed 'local' deviations to preserve active set
% Error between analytical and numerical (actual/true) solutions

clear;
clc;

%% Nominal CPT parameter values 

% Within parameter bounds
alpha = 0.82;
beta = 0.8;
lambda = 2.25;
p = 0.75; % Also try w/ p = 0.25
nominal = [alpha beta lambda p];

%%
num_devs = 15;
devs = linspace(-20,20,num_devs); 
alpha_vals = alpha + alpha.*(devs./100);
beta_vals = beta + beta.*(devs./100);
lambda_vals = lambda + lambda.*(devs./100);
p_vals = p + p.*(devs./100);

%% Loop through several possible scenarios
% Vary u0, b_sm, lb and ub in discrete manner

u0_lb = -10;
u0_ub = 10;
b_sm_lb = -0.99;
b_sm_ub = -0.01;
p_lb1 = 1;
p_lb2 = 5;
p_ub1 = 6;
p_ub2 = 20;

num_scenarios = 20;
result.u0 = zeros(1,num_scenarios);
result.b_sm = zeros(1,num_scenarios);
result.x1 = zeros(1,num_scenarios);
result.x2 = zeros(1,num_scenarios);
result.valid = zeros(1,num_scenarios);

% tariff bounds
result.lb = zeros(1,num_scenarios); 
result.ub = zeros(1,num_scenarios);

% error between analytical and numerical solutions
result.error_soln = zeros(num_scenarios,4);
result.error_cost = zeros(num_scenarios,4);
result.error_mismatch = zeros(num_scenarios,4);
result.max_devs = zeros(num_scenarios,4);
result.dGamma_dTheta = zeros(num_scenarios,4);
result.dF_dTheta = zeros(num_scenarios,4);

for i = 1:num_scenarios 
    u0 = u0_lb + (u0_ub - u0_lb)*rand;
    b_sm = b_sm_lb + (b_sm_ub - b_sm_lb)*rand;
    lb = p_lb1 + (p_lb2 - p_lb1)*rand;
    ub = p_ub1 + (p_ub2 - p_ub1)*rand;
    
    % Determine limiting values for x1 & x2 that ensure non-trivial choice scenario
    x1_ub = u0 - b_sm*lb;
    x2_lb = u0 - b_sm*ub;

    % Then randomly pick values
    x1 = x1_ub - abs(x1_ub)*rand; 
    x2 = x2_lb + abs(x2_lb)*rand;

    % Sanity check
    u1_max = x1 + b_sm*lb;
    u2_min = x2 + b_sm*ub;
    result.valid(i) = (u1_max <= u0 && u0 <= u2_min);
  
    scenario_vals.u0 = u0;
    scenario_vals.x1 = x1;
    scenario_vals.x2 = x2;
    scenario_vals.b_sm = b_sm;
    scenario_vals.R = 'Best case';
    scenario_vals.lb = lb;
    scenario_vals.ub = ub;

    result.u0(i) = u0;
    result.b_sm(i) = b_sm;
    result.lb(i) = lb;
    result.ub(i) = ub; 
    result.x1(i) = x1;
    result.x2(i) = x2;
    params = {'alpha','beta','lambda','p'};
    
    for j = 1:4
        theta = params{j};
        [opt_soln,opt_cost,~,~,mismatch_loss,~,~] = derivations_numerical(theta,nominal,scenario_vals,devs,0);
        [dGamma_dTheta_nom,dF_dTheta_nom,opt_soln_an,opt_cost_an,mismatch_loss_an,~,max_dev,~,~] = derivations_local(theta,nominal,scenario_vals,devs,0);
        str = join([theta,num2str(i)]);
        str_sol = join([str,'opt_soln']);
        str_sol_an = join([str,'opt_soln_an']);
        str_cost = join([str,'opt_cost']);
        str_cost_an = join([str,'opt_cost_an']);
        str_mismatch = join([str,'mismatch']);
        str_mismatch_an = join([str,'mismatch_an']);
        result.(str_sol) = opt_soln;
        result.(str_cost) = opt_cost;
        result.(str_mismatch) = mismatch_loss;
        result.(str_sol_an) = opt_soln_an;
        result.(str_cost_an) = opt_cost_an;
        result.(str_mismatch_an) = mismatch_loss_an;
        
        % Relative (%) errors
        result.error_soln(i,j) = (norm(opt_soln - opt_soln_an)/norm(opt_soln))*100;
        result.error_cost(i,j) = (norm(opt_cost - opt_cost_an)/norm(opt_cost))*100;
        result.error_mismatch(i,j) = (norm(mismatch_loss - mismatch_loss_an)/norm(mismatch_loss))*100;
    
        result.max_devs(i,j) = max_dev;
        result.dGamma_dTheta(i,j) = dGamma_dTheta_nom;
        result.dF_dTheta(i,j) = dF_dTheta_nom;
    end
end
%% Make and save plots
load('C:\Users\vinee\Dropbox\CDC 2020_Vineet and Yue\CDC_2020_Results\Test_20RandomScenarios_v4_p0.25.mat')
num_scenarios = 20;

for i = 1:num_scenarios
    str_sol = join(['OptSoln',num2str(i)]);
    str_cost = join(['OptCost',num2str(i)]);
    str_mismatch = join(['Mismatch',num2str(i)]);
    str_compare_sol = join(['SolComparison',num2str(i)]);
    str_compare_cost = join(['CostComparison',num2str(i)]);
    directory = 'C:\Users\vinee\Dropbox\CDC 2020_Vineet and Yue\CDC_2020_Results\Test_20RandomScenarios_v4_p0.25';
    
    str_alpha = join(['alpha', num2str(i)]);
    str_sol_alpha = join([str_alpha,'opt_soln']);
    str_cost_alpha = join([str_alpha,'opt_cost']);
    str_mismatch_alpha = join([str_alpha,'mismatch']);
    str_sol_alpha_an = join([str_alpha,'opt_soln_an']);
    str_cost_alpha_an = join([str_alpha,'opt_cost_an']);
    str_mismatch_alpha_an = join([str_alpha,'mismatch_an']);    
    
    str_beta = join(['beta', num2str(i)]);
    str_sol_beta = join([str_beta,'opt_soln']);
    str_cost_beta = join([str_beta,'opt_cost']);
    str_mismatch_beta = join([str_beta,'mismatch']);
    str_sol_beta_an = join([str_beta,'opt_soln_an']);
    str_cost_beta_an = join([str_beta,'opt_cost_an']);
    str_mismatch_beta_an = join([str_beta,'mismatch_an']);
    
    str_lambda = join(['lambda', num2str(i)]);
    str_sol_lambda = join([str_lambda,'opt_soln']);
    str_cost_lambda = join([str_lambda,'opt_cost']);
    str_mismatch_lambda = join([str_lambda,'mismatch']);
    str_sol_lambda_an = join([str_lambda,'opt_soln_an']);
    str_cost_lambda_an = join([str_lambda,'opt_cost_an']);
    str_mismatch_lambda_an = join([str_lambda,'mismatch_an']);
    
    str_p = join(['p', num2str(i)]);
    str_sol_p = join([str_p,'opt_soln']);
    str_cost_p = join([str_p,'opt_cost']);
    str_mismatch_p = join([str_p,'mismatch']);
    str_sol_p_an = join([str_p,'opt_soln_an']);
    str_cost_p_an = join([str_p,'opt_cost_an']);
    str_mismatch_p_an = join([str_p,'mismatch_an']);
    
    % Optimal solution plots - numerical
    ax1(1) = subplot(2,2,1);
    hold on
    plot(alpha_vals,result.(str_sol_alpha));
    plot(alpha_vals(8)*ones(1,15),linspace(0,result.(str_sol_alpha)(8),15),'--')
    plot(alpha_vals(8),result.(str_sol_alpha)(8),'x','LineWidth',1.25)
    xlabel('\alpha');
    ylabel('Optimal solution ($)');

    ax1(2) = subplot(2,2,2);
    plot(beta_vals,result.(str_sol_beta));
    xlabel('\beta');
    ylabel('Optimal solution ($)');

    ax1(3) = subplot(2,2,3);
    plot(lambda_vals,result.(str_sol_lambda));
    xlabel('\lambda');
    ylabel('Optimal solution ($)');
    
    ax1(4) = subplot(2,2,4);
    plot(p_vals,result.(str_sol_p));
    xlabel('p');
    ylabel('Optimal solution ($)');
   
    linkaxes(ax1,'y');
    ax1(1).YLim = [result.lb(i) result.ub(i)];
    
    sgtitle(sprintf('u_0=%0.2f, b_{sm}=%0.2f, lb=%0.2f, ub=%0.2f, x_1=%0.2f, x_2=%0.2f',...
        result.u0(i),result.b_sm(i),result.lb(i),result.ub(i),result.x1(i),result.x2(i)))
    name = fullfile(directory,str_sol);
    saveas(gcf,name,'png')
    close
    
    % Optimal objective/value function plots - numerical
    ax2(1) = subplot(2,2,1);
    plot(alpha_vals,result.(str_cost_alpha));
    xlabel('\alpha');
    ylabel('Optimal objective ($)');

    ax2(2) = subplot(2,2,2);
    plot(beta_vals,result.(str_cost_beta));
    xlabel('\beta');
    ylabel('Optimal objective ($)');

    ax2(3) = subplot(2,2,3);
    plot(lambda_vals,result.(str_cost_lambda));
    xlabel('\lambda');
    ylabel('Optimal objective ($)');
    
    ax2(4) = subplot(2,2,4);
    plot(p_vals,result.(str_cost_p));
    xlabel('p');
    ylabel('Optimal objective ($)');
    
    linkaxes(ax2,'y');
    sgtitle(sprintf('u_0=%0.2f, b_{sm}=%0.2f, lb=%0.2f, ub=%0.2f, x_1=%0.2f, x_2=%0.2f',...
        result.u0(i),result.b_sm(i),result.lb(i),result.ub(i),result.x1(i),result.x2(i)))
    name = fullfile(directory,str_cost);
    saveas(gcf,name,'png')
    close
    
    % Mismatch loss plots - numerical
    ax3(1) = subplot(2,2,1);
    plot(alpha_vals,result.(str_mismatch_alpha));
    xlabel('\alpha');
    ylabel('Mismatch loss ($)');

    ax3(2) = subplot(2,2,2);
    plot(beta_vals,result.(str_mismatch_beta));
    xlabel('\beta');
    ylabel('Mismatch loss ($)');

    ax3(3) = subplot(2,2,3);
    plot(lambda_vals,result.(str_mismatch_lambda));
    xlabel('\lambda');
    ylabel('Mismatch loss ($)');
    
    ax3(4) = subplot(2,2,4);
    plot(p_vals,result.(str_mismatch_p));
    xlabel('p');
    ylabel('Mismatch loss ($)');
    
    linkaxes(ax3,'y');
    sgtitle(sprintf('u_0=%0.2f, b_{sm}=%0.2f, lb=%0.2f, ub=%0.2f, x_1=%0.2f, x_2=%0.2f',...
        result.u0(i),result.b_sm(i),result.lb(i),result.ub(i),result.x1(i),result.x2(i)))
    name = fullfile(directory,str_mismatch);
    saveas(gcf,name,'png')
    close

    % Optimal solution plots - comparison
    ax4(1) = subplot(2,2,1);
    hold on
    plot(alpha_vals,result.(str_sol_alpha));
    plot(alpha_vals,result.(str_sol_alpha_an));
    legend('Numerical','Analytical','Location','best')
    xlabel('\alpha');
    ylabel('Optimal solution ($)');
    hold off

    ax4(2) = subplot(2,2,2);
    hold on
    plot(beta_vals,result.(str_sol_beta));
    plot(beta_vals,result.(str_sol_beta_an));
    legend('Numerical','Analytical','Location','best')
    xlabel('\beta');
    ylabel('Optimal solution ($)');
    hold off

    ax4(3) = subplot(2,2,3);
    hold on
    plot(lambda_vals,result.(str_sol_lambda));
    plot(lambda_vals,result.(str_sol_lambda_an));
    legend('Numerical','Analytical','Location','best')
    xlabel('\lambda');
    ylabel('Optimal solution ($)');
    hold off
    
    ax4(4) = subplot(2,2,4);
    hold on
    plot(p_vals,result.(str_sol_p));
    plot(p_vals,result.(str_sol_p_an));
    legend('Numerical','Analytical','Location','best')
    xlabel('p');
    ylabel('Optimal solution ($)');
    hold off
   
    linkaxes(ax4,'y');
    ax4(1).YLim = [result.lb(i) result.ub(i)];
    
    sgtitle(sprintf('u_0=%0.2f, b_{sm}=%0.2f, lb=%0.2f, ub=%0.2f, x_1=%0.2f, x_2=%0.2f',...
        result.u0(i),result.b_sm(i),result.lb(i),result.ub(i),result.x1(i),result.x2(i)))
    name = fullfile(directory,str_compare_sol);
    saveas(gcf,name,'png')
    close
    
    % Optimal objective plots - comparison
    ax5(1) = subplot(2,2,1);
    hold on
    plot(alpha_vals,result.(str_cost_alpha));
    plot(alpha_vals,result.(str_cost_alpha_an));
    legend('Numerical','Analytical','Location','best')
    xlabel('\alpha');
    ylabel('Optimal objective ($)');
    hold off

    ax5(2) = subplot(2,2,2);
    hold on
    plot(beta_vals,result.(str_cost_beta));
    plot(beta_vals,result.(str_cost_beta_an));
    legend('Numerical','Analytical','Location','best')
    xlabel('\beta');
    ylabel('Optimal objective ($)');
    hold off

    ax5(3) = subplot(2,2,3);
    hold on
    plot(lambda_vals,result.(str_cost_lambda));
    plot(lambda_vals,result.(str_cost_lambda_an));
    legend('Numerical','Analytical','Location','best')
    xlabel('\lambda');
    ylabel('Optimal objective ($)');
    hold off
    
    ax5(4) = subplot(2,2,4);
    hold on
    plot(p_vals,result.(str_cost_p));
    plot(p_vals,result.(str_cost_p_an));
    legend('Numerical','Analytical','Location','best')
    xlabel('p');
    ylabel('Optimal objective ($)');
    hold off
   
    linkaxes(ax5,'y');    
    sgtitle(sprintf('u_0=%0.2f, b_{sm}=%0.2f, lb=%0.2f, ub=%0.2f, x_1=%0.2f, x_2=%0.2f',...
        result.u0(i),result.b_sm(i),result.lb(i),result.ub(i),result.x1(i),result.x2(i)))
    name = fullfile(directory,str_compare_cost);
    saveas(gcf,name,'png')
    close
end

