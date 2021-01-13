load('C:\Users\vinee\Dropbox\CDC 2020_Vineet and Yue\CDC_2020_Results\Test_20RandomScenarios_v4_trial3')
num_scenarios = 20;

for i = 12%:num_scenarios
    str_sol = join(['OptSoln',num2str(i)]);
    str_cost = join(['OptCost',num2str(i)]);
    str_mismatch = join(['Mismatch',num2str(i)]);
    str_compare_sol = join(['SolComparison',num2str(i)]);
    str_compare_cost = join(['CostComparison',num2str(i)]);
    directory = 'C:\Users\vinee\Dropbox\CDC 2020_Vineet and Yue\CDC_2020_Results\Final_plots_v3';
    
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
    set(gca,'FontSize',12)
    hold on
    plot(alpha_vals,result.(str_sol_alpha),'LineWidth',1.5);
    plot(alpha_vals(8)*ones(1,15),linspace(0,result.(str_sol_alpha)(8),15),'--','LineWidth',1.5)
    plot(alpha_vals(8),result.(str_sol_alpha)(8),'x','LineWidth',1.5)
    xlabel('\alpha');
    ylabel('Optimal solution ($)');
    xlim([alpha_vals(1) alpha_vals(end)]);
    %ylim([7 15])

    ax1(2) = subplot(2,2,2);
    set(gca,'FontSize',12)
    hold on
    plot(beta_vals,result.(str_sol_beta),'LineWidth',1.5);
    plot(beta_vals(8)*ones(1,15),linspace(0,result.(str_sol_beta)(8),15),'--','LineWidth',1.5)
    plot(beta_vals(8),result.(str_sol_beta)(8),'x','LineWidth',1.5)
    xlabel('\beta');
    ylabel('Optimal solution ($)');
    xlim([beta_vals(1) beta_vals(end)]);
    %ylim([7 15])

    ax1(3) = subplot(2,2,3);
    set(gca,'FontSize',12)
    hold on
    plot(lambda_vals,result.(str_sol_lambda),'LineWidth',1.5);
    plot(lambda_vals(8)*ones(1,15),linspace(0,result.(str_sol_lambda)(8),15),'--','LineWidth',1.5)
    plot(lambda_vals(8),result.(str_sol_lambda)(8),'x','LineWidth',1.5)
    xlabel('\lambda');
    ylabel('Optimal solution ($)');
    xlim([lambda_vals(1) lambda_vals(end)]);
    %ylim([7 15])
    
    ax1(4) = subplot(2,2,4);
    set(gca,'FontSize',12)
    hold on
    plot(p_vals,result.(str_sol_p),'LineWidth',1.5);
    plot(p_vals(8)*ones(1,15),linspace(0,result.(str_sol_p)(8),15),'--','LineWidth',1.5)
    plot(p_vals(8),result.(str_sol_p)(8),'x','LineWidth',1.5)
    xlabel('$p$','Interpreter','latex');
    ylabel('Optimal solution ($)');
    xlim([p_vals(1) p_vals(end)]);
    %ylim([7 15])
    
    linkaxes(ax1,'y');    
    name = fullfile(directory,str_sol);
    saveas(gcf,name,'png')
    close
    
    % Optimal objective/value function plots - numerical
    ax2(1) = subplot(2,2,1);
    set(gca,'FontSize',12)
    hold on
    plot(alpha_vals,result.(str_cost_alpha),'LineWidth',1.5);
    plot(alpha_vals(8)*ones(1,15),linspace(0,result.(str_cost_alpha)(8),15),'--','LineWidth',1.5)
    plot(alpha_vals(8),result.(str_cost_alpha)(8),'x','LineWidth',1.5)
    xlabel('\alpha');
    ylabel('Optimal objective ($)');
    xlim([alpha_vals(1) alpha_vals(end)]);
    %ylim([1 6])
    
    ax2(2) = subplot(2,2,2);
    set(gca,'FontSize',12)
    hold
    plot(beta_vals,result.(str_cost_beta),'LineWidth',1.5);
    plot(beta_vals(8)*ones(1,15),linspace(0,result.(str_cost_beta)(8),15),'--','LineWidth',1.5)
    plot(beta_vals(8),result.(str_cost_beta)(8),'x','LineWidth',1.5)
    xlabel('\beta');
    ylabel('Optimal objective ($)');
    xlim([beta_vals(1) beta_vals(end)]);
    %ylim([1 6])

    ax2(3) = subplot(2,2,3);
    set(gca,'FontSize',12)
    hold on
    plot(lambda_vals,result.(str_cost_lambda),'LineWidth',1.5);
    plot(lambda_vals(8)*ones(1,15),linspace(0,result.(str_cost_lambda)(8),15),'--','LineWidth',1.5)
    plot(lambda_vals(8),result.(str_cost_lambda)(8),'x','LineWidth',1.5)
    xlabel('\lambda');
    ylabel('Optimal objective ($)');
    xlim([lambda_vals(1) lambda_vals(end)]);
    %ylim([1 6])
    
    ax2(4) = subplot(2,2,4);
    set(gca,'FontSize',12)
    hold on
    plot(p_vals,result.(str_cost_p),'LineWidth',1.5);
    plot(p_vals(8)*ones(1,15),linspace(0,result.(str_cost_p)(8),15),'--','LineWidth',1.5)
    plot(p_vals(8),result.(str_cost_p)(8),'x','LineWidth',1.5)
    xlabel('$p$','Interpreter','latex');
    ylabel('Optimal objective ($)');
    xlim([p_vals(1) p_vals(end)]);
    %ylim([1 6])
    
    linkaxes(ax2,'y');
    name = fullfile(directory,str_cost);
    saveas(gcf,name,'png')
    close
    
    % Mismatch loss plots - numerical
    ax3(1) = subplot(2,2,1);
    set(gca,'FontSize',12)
    hold on
    plot(alpha_vals,result.(str_mismatch_alpha),'LineWidth',1.5);
    plot(alpha_vals(8)*ones(1,15),linspace(0,result.(str_mismatch_alpha)(8),15),'--','LineWidth',1.5)
    plot(alpha_vals(8),result.(str_mismatch_alpha)(8),'x','LineWidth',1.5)
    xlabel('\alpha');
    ylabel('Mismatch loss ($)');
    xlim([alpha_vals(1) alpha_vals(end)]);

    ax3(2) = subplot(2,2,2);
    set(gca,'FontSize',12)
    hold on
    plot(beta_vals,result.(str_mismatch_beta),'LineWidth',1.5);
    plot(beta_vals(8)*ones(1,15),linspace(0,result.(str_mismatch_beta)(8),15),'--','LineWidth',1.5)
    plot(beta_vals(8),result.(str_mismatch_beta)(8),'x','LineWidth',1.5)
    xlabel('\beta');
    ylabel('Mismatch loss ($)');
    xlim([beta_vals(1) beta_vals(end)]);

    ax3(3) = subplot(2,2,3);
    set(gca,'FontSize',12)
    hold on
    plot(lambda_vals,result.(str_mismatch_lambda),'LineWidth',1.5);
    plot(lambda_vals(8)*ones(1,15),linspace(0,result.(str_mismatch_lambda)(8),15),'--','LineWidth',1.5)
    plot(lambda_vals(8),result.(str_mismatch_lambda)(8),'x','LineWidth',1.5)
    xlabel('\lambda');
    ylabel('Mismatch loss ($)');
    xlim([lambda_vals(1) lambda_vals(end)]);
    
    ax3(4) = subplot(2,2,4);
    set(gca,'FontSize',12)
    hold on
    plot(p_vals,result.(str_mismatch_p),'LineWidth',1.5);
    plot(p_vals(8)*ones(1,15),linspace(0,result.(str_mismatch_p)(8),15),'--','LineWidth',1.5)
    plot(p_vals(8),result.(str_mismatch_p)(8),'x','LineWidth',1.5)
    xlabel('$p$','Interpreter','latex');
    ylabel('Mismatch loss ($)');
    xlim([p_vals(1) p_vals(end)]);
    
    linkaxes(ax3,'y');
    name = fullfile(directory,str_mismatch);
    saveas(gcf,name,'png')
    close

    % Optimal solution plots - comparison
    ax4(1) = subplot(2,2,1);
    set(gca,'FontSize',12)
    hold on
    plot(alpha_vals,result.(str_sol_alpha),'LineWidth',1.25);
    plot(alpha_vals,result.(str_sol_alpha_an),'LineWidth',1.25)
    legend('Numerical','Analytical','Location','best','AutoUpdate','off')
    plot(alpha_vals(8)*ones(1,15),linspace(0,result.(str_sol_alpha)(8),15),'--','LineWidth',1.25)
    plot(alpha_vals(8),result.(str_sol_alpha)(8),'x','LineWidth',1.25)
    xlabel('\alpha');
    ylabel('Optimal solution ($)');
    hold off
    xlim([alpha_vals(1) alpha_vals(end)]);

    ax4(2) = subplot(2,2,2);
    set(gca,'FontSize',12)
    hold on
    plot(beta_vals,result.(str_sol_beta),'LineWidth',1.25);
    plot(beta_vals,result.(str_sol_beta_an),'LineWidth',1.25);
    legend('Numerical','Analytical','Location','best','AutoUpdate','off')
    plot(beta_vals(8)*ones(1,15),linspace(0,result.(str_sol_beta)(8),15),'--','LineWidth',1.25)
    plot(beta_vals(8),result.(str_sol_beta)(8),'x','LineWidth',1.25)
    xlabel('\beta');
    ylabel('Optimal solution ($)');
    hold off
    xlim([beta_vals(1) beta_vals(end)]);

    ax4(3) = subplot(2,2,3);
    set(gca,'FontSize',12)
    hold on
    plot(lambda_vals,result.(str_sol_lambda),'LineWidth',1.25);
    plot(lambda_vals,result.(str_sol_lambda_an),'LineWidth',1.25);
    legend('Numerical','Analytical','Location','best','AutoUpdate','off')
    plot(lambda_vals(8)*ones(1,15),linspace(0,result.(str_sol_lambda)(8),15),'--','LineWidth',1.25)
    plot(lambda_vals(8),result.(str_sol_lambda)(8),'x','LineWidth',1.25)
    xlabel('\lambda');
    ylabel('Optimal solution ($)');
    hold off
    xlim([lambda_vals(1) lambda_vals(end)]);
    
    ax4(4) = subplot(2,2,4);
    set(gca,'FontSize',12)
    hold on
    plot(p_vals,result.(str_sol_p),'LineWidth',1.25);
    plot(p_vals,result.(str_sol_p_an),'LineWidth',1.25);
    legend('Numerical','Analytical','Location','best','AutoUpdate','off')
    plot(p_vals(8)*ones(1,15),linspace(0,result.(str_sol_p)(8),15),'--','LineWidth',1.25)
    plot(p_vals(8),result.(str_sol_p)(8),'x','LineWidth',1.25)
    xlabel('$p$','Interpreter','latex');
    ylabel('Optimal solution ($)');
    hold off
    xlim([p_vals(1) p_vals(end)]);
   
    linkaxes(ax4,'y');  
    name = fullfile(directory,str_compare_sol);
    saveas(gcf,name,'png')
    close
    
    % Optimal objective plots - comparison
    ax5(1) = subplot(2,2,1);
    set(gca,'FontSize',12)
    hold on
    plot(alpha_vals,result.(str_cost_alpha),'LineWidth',1.25);
    plot(alpha_vals,result.(str_cost_alpha_an),'LineWidth',1.25);
    legend('Numerical','Analytical','Location','best','AutoUpdate','off')
    plot(alpha_vals(8)*ones(1,15),linspace(0,result.(str_cost_alpha)(8),15),'--','LineWidth',1.25)
    plot(alpha_vals(8),result.(str_cost_alpha)(8),'x','LineWidth',1.25)
    xlabel('\alpha');
    ylabel('Optimal objective ($)');
    hold off
    xlim([alpha_vals(1) alpha_vals(end)]);
    
    ax5(2) = subplot(2,2,2);
    set(gca,'FontSize',12)
    hold on
    plot(beta_vals,result.(str_cost_beta),'LineWidth',1.25);
    plot(beta_vals,result.(str_cost_beta_an),'LineWidth',1.25);
    legend('Numerical','Analytical','Location','best','AutoUpdate','off')
    plot(beta_vals(8)*ones(1,15),linspace(0,result.(str_cost_beta)(8),15),'--','LineWidth',1.25)
    plot(beta_vals(8),result.(str_cost_beta)(8),'x','LineWidth',1.25)
    xlabel('\beta');
    ylabel('Optimal objective ($)');
    hold off
    xlim([beta_vals(1) beta_vals(end)]);

    ax5(3) = subplot(2,2,3);
    set(gca,'FontSize',12)
    hold on
    plot(lambda_vals,result.(str_cost_lambda),'LineWidth',1.25);
    plot(lambda_vals,result.(str_cost_lambda_an),'LineWidth',1.25);
    legend('Numerical','Analytical','Location','best','AutoUpdate','off')
    plot(lambda_vals(8)*ones(1,15),linspace(0,result.(str_cost_lambda)(8),15),'--','LineWidth',1.25)
    plot(lambda_vals(8),result.(str_cost_lambda)(8),'x','LineWidth',1.25)
    xlabel('\lambda');
    ylabel('Optimal objective ($)');
    hold off
    xlim([lambda_vals(1) lambda_vals(end)]);
    
    ax5(4) = subplot(2,2,4);
    set(gca,'FontSize',12)
    hold on
    plot(p_vals,result.(str_cost_p),'LineWidth',1.25);
    plot(p_vals,result.(str_cost_p_an),'LineWidth',1.25);
    legend('Numerical','Analytical','Location','best','AutoUpdate','off')
    plot(p_vals(8)*ones(1,15),linspace(0,result.(str_cost_p)(8),15),'--','LineWidth',1.25)
    plot(p_vals(8),result.(str_cost_p)(8),'x','LineWidth',1.25)
    xlabel('$p$','Interpreter','latex');
    ylabel('Optimal objective ($)');
    hold off
    xlim([p_vals(1) p_vals(end)]);
   
    linkaxes(ax5,'y');    
    name = fullfile(directory,str_compare_cost);
    saveas(gcf,name,'png')
    close
end