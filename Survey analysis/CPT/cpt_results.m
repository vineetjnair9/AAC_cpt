%% Loading data
T = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Data\Panel\Full_launch_6-10_n955');
table_orig = T.Fulllaunch610n955;
table = convertvars(table_orig,[7:8,10:12,48:139],'double');

load('valid_indices.mat')

% Select subset of table rows corresponding to valid respondents
table = table(valid_indices,:);

num = length(valid_indices);
respondents = 1:1:num;
respondents = respondents';

[I,I_id] = findgroups(table.income);

%% Convert z_final array to table
z = array2table(z_final(:,3:10),'VariableNames',...
    {'Respondent','nu','alpha_plus','alpha_minus','beta_plus','beta_minus','lambda','error'});

clear T z_final table_orig

%% Scatter plots
reg_vals = 1.5;
for i = 1:length(reg_vals)
    reg = reg_vals(i);
    rows = (round(z.nu,1) == reg);
    cpt = table2array(z(rows,3:7));
    error = table2array(z(rows,8));

    % cpt = cpt_lottery;
    alpha_plus = cpt(:,1);
    alpha_minus = cpt(:,2);
    beta_plus = cpt(:,3);
    beta_minus = cpt(:,4);
    lambda = cpt(:,5).*100;

    % Scatter plots
    dir = 'C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\SM Thesis\Figures\CPT_panel_reg_trials\';
    weight = 'OG_';
    norm = 'trueCE';

    figure(1)
    scatter(respondents,alpha_plus)
    ylabel('$\alpha^+$','Interpreter','latex')
    xlabel('Respondent','Interpreter','latex')
    name = join([dir,'cpt_panel_final_','dynSMoDS_',weight,norm,'_reg_',num2str(reg),'_scatter_alpha+.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close

    figure(2)
    scatter(respondents,alpha_minus)
    ylabel('$\alpha^-$','Interpreter','latex')
    xlabel('Respondent','Interpreter','latex')
    name = join([dir,'cpt_panel_final_','dynSMoDS_',weight,norm,'_reg_',num2str(reg),'_scatter_alpha-.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close

    figure(3)
    scatter(respondents,beta_plus)
    ylabel('$\beta^+$','Interpreter','latex')
    xlabel('Respondent','Interpreter','latex')
    name = join([dir,'cpt_panel_final_','dynSMoDS_',weight,norm,'_reg_',num2str(reg),'_scatter_beta+.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close

    figure(4)
    scatter(respondents,beta_minus)
    ylabel('$\beta^-$','Interpreter','latex')
    xlabel('Respondent','Interpreter','latex')
    name = join([dir,'cpt_panel_final_','dynSMoDS_',weight,norm,'_reg_',num2str(reg),'_scatter_beta-.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close

    figure(5)
    scatter(respondents,lambda)
    ylabel('$\lambda$','Interpreter','latex')
    xlabel('Respondent','Interpreter','latex')
    name = join([dir,'cpt_panel_final_','dynSMoDS_',weight,norm,'_reg_',num2str(reg),'_scatter_lambda.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close
    
    figure(6)
    scatter(respondents,error)
    ylabel('Squared norm of residual','Interpreter','latex')
    xlabel('Respondent','Interpreter','latex')
    name = join([dir,'cpt_panel_final_',weight,norm,'_reg_',num2str(reg),'_scatter_error.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close
end

%% Histograms
reg_vals = 0:0.5:3;
% Scatter plots
dir = 'C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\SM Thesis\Figures\CPT_panel_reg_trials\';
weight = 'OG_';
norm = 'trueCE';
R = 'dynSMoDS_';
for i = 1:length(reg_vals)
    reg = reg_vals(i);
    rows = (round(z.nu,1) == reg);
    cpt = table2array(z(rows,3:7));
    error = table2array(z(rows,8));

    % cpt = cpt_lottery;
    alpha_plus = cpt(:,1);
    alpha_minus = cpt(:,2);
    beta_plus = cpt(:,3);
    beta_minus = cpt(:,4);
    lambda = cpt(:,5).*100;

    figure(1)
    histogram(alpha_plus)
    xlabel('$\alpha^+$','Interpreter','latex')
    ylabel('Number of respondents','Interpreter','latex')
    name = join([dir,'cpt_panel_final_',R,weight,norm,'_reg_',num2str(reg),'_hist_alpha+.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close

    figure(2)
    histogram(alpha_minus)
    xlabel('$\alpha^-$','Interpreter','latex')
    ylabel('Number of respondents','Interpreter','latex')
    name = join([dir,'cpt_panel_final_',R,weight,norm,'_reg_',num2str(reg),'_hist_alpha-.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close

    figure(3)
    histogram(beta_plus)
    xlabel('$\beta^+$','Interpreter','latex')
    ylabel('Number of respondents','Interpreter','latex')
    name = join([dir,'cpt_panel_final_',R,weight,norm,'_reg_',num2str(reg),'_hist_beta+.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close

    figure(4)
    histogram(beta_minus)
    xlabel('$\beta^-$','Interpreter','latex')
    ylabel('Number of respondents','Interpreter','latex')
    name = join([dir,'cpt_panel_final_',R,weight,norm,'_reg_',num2str(reg),'_hist_beta-.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close

    figure(5)
    histogram(lambda)
    xlabel('$\lambda$','Interpreter','latex')
    ylabel('Number of respondents','Interpreter','latex')
    name = join([dir,'cpt_panel_final_',R,weight,norm,'_reg_',num2str(reg),'_hist_lambda.png']);
    exportgraphics(gcf,name,'Resolution',300)
    close
end

%%
regs = 1.5;
rows = (round(z.nu,1) == reg);
cpt = table2array(z(rows,3:7));
error = table2array(z(rows,8));

lb = [0,0,0,0,0.01];
ub = [1,1,1,1,1];
% cpt = z_final(:,5:9);
lb_error_norms = vecnorm(cpt - lb,2,2).^2;
ub_error_norms = vecnorm(cpt - ub,2,2).^2;
% regs = z_final(:,4);
adjusted_errors = error - regs.*lb_error_norms - regs.*ub_error_norms;
%%
cpt = cpt_lottery;
alpha_plus = cpt(:,1);
alpha_minus = cpt(:,2);
beta_plus = cpt(:,3);
beta_minus = cpt(:,4);
lambda = cpt(:,5);

%% Population averages
% Mean
alpha_plus_mean = mean(alpha_plus)
alpha_minus_mean = mean(alpha_minus)
beta_plus_mean = mean(beta_plus)
beta_minus_mean = mean(beta_minus)
lambda_mean = mean(lambda)
error_mean = mean(error)

%% Median
clc;
alpha_plus_median = median(alpha_plus)
alpha_minus_median = median(alpha_minus)
beta_plus_median = median(beta_plus)
beta_minus_median = median(beta_minus)
lambda_median = median(lambda)
error_median = median(error)

%% Minimum values
clc;
alpha_plus_min = min(alpha_plus)
alpha_minus_min = min(alpha_minus)
beta_plus_min = min(beta_plus)
beta_minus_min = min(beta_minus)
lambda_min = min(lambda)

%% Maximum values
clc;
alpha_plus_max = max(alpha_plus)
alpha_minus_max = max(alpha_minus)
beta_plus_max = max(beta_plus)
beta_minus_max = max(beta_minus)
lambda_max = max(lambda)

%% Standard deviations
clc;
alpha_plus_sd = std(alpha_plus)
alpha_minus_sd = std(alpha_minus)
beta_plus_sd = std(beta_plus)
beta_minus_sd = std(beta_minus)
lambda_sd = std(lambda)
error_sd = std(error)

%% Standard errors
clc;
alpha_plus_sd = std(alpha_plus)/sqrt(955)
alpha_minus_sd = std(alpha_minus)/sqrt(955)
beta_plus_sd = std(beta_plus)/sqrt(955)
beta_minus_sd = std(beta_minus)/sqrt(955)
lambda_sd = std(lambda)/sqrt(955)

%% Box plots 
dir = 'C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\SM Thesis\Figures\CPT_demographic\';
demo = 'gender';
% Gender
figure(1)
boxplot(alpha_plus,table.gender,'labels',{'Male','Female','Non-binary'})
ylabel('\alpha^+')
name = join([dir,demo,'alpha+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(2)
boxplot(alpha_minus,table.gender,'labels',{'Male','Female','Non-binary'})
ylabel('\alpha^-')
name = join([dir,demo,'alpha-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(3)
boxplot(beta_plus,table.gender,'labels',{'Male','Female','Non-binary'})
ylabel('\beta^+')
name = join([dir,demo,'beta+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(4)
boxplot(beta_minus,table.gender,'labels',{'Male','Female','Non-binary'})
ylabel('\beta^-')
name = join([dir,demo,'beta-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(5)
boxplot(lambda,table.gender,'labels',{'Male','Female','Non-binary'})
ylabel('\lambda')
name = join([dir,demo,'lambda.png']);
exportgraphics(gcf,name,'Resolution',300)
close

%% Income
demo = 'income';
Income_order = string([I_id(11); I_id(1); I_id(3:10); I_id(2); I_id(12)]);
figure(1)
boxplot(alpha_plus,table.income,'PlotStyle','compact','GroupOrder',Income_order,...
    'labels',{'< $10K','$10-20K','$20-30K','$30-40K','$40-50K','$50-60K','$60-70K','$70-80K','$80-90K','$90-100K','$100-150K','>$150K'})
ylabel('\alpha^+')
name = join([dir,demo,'alpha+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(2)
boxplot(alpha_minus,table.income,'PlotStyle','compact','GroupOrder',Income_order,...
        'labels',{'< $10K','$10-20K','$20-30K','$30-40K','$40-50K','$50-60K','$60-70K','$70-80K','$80-90K','$90-100K','$100-150K','>$150K'})
ylabel('\alpha^-')
name = join([dir,demo,'alpha-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(3)
boxplot(beta_plus,table.income,'PlotStyle','compact','GroupOrder',Income_order,...
        'labels',{'< $10K','$10-20K','$20-30K','$30-40K','$40-50K','$50-60K','$60-70K','$70-80K','$80-90K','$90-100K','$100-150K','>$150K'})
ylabel('\beta^+')
name = join([dir,demo,'beta+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(4)
boxplot(beta_minus,table.income,'PlotStyle','compact','GroupOrder',Income_order,...
        'labels',{'< $10K','$10-20K','$20-30K','$30-40K','$40-50K','$50-60K','$60-70K','$70-80K','$80-90K','$90-100K','$100-150K','>$150K'})
ylabel('\beta^-')
name = join([dir,demo,'beta-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(5)
boxplot(lambda,table.income,'PlotStyle','compact','GroupOrder',Income_order,...
        'labels',{'< $10K','$10-20K','$20-30K','$30-40K','$40-50K','$50-60K','$60-70K','$70-80K','$80-90K','$90-100K','$100-150K','>$150K'})
ylabel('\lambda')
name = join([dir,demo,'lambda.png']);
exportgraphics(gcf,name,'Resolution',300)
close

%% Age
demo = 'age';
figure(1)
boxplot(alpha_plus,table.age)
ylabel('\alpha^+')
name = join([dir,demo,'alpha+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(2)
boxplot(alpha_minus,table.age)
ylabel('\alpha^-')
name = join([dir,demo,'alpha-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(3)
boxplot(beta_plus,table.age)
ylabel('\beta^+')
name = join([dir,demo,'beta+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(4)
boxplot(beta_minus,table.age)
ylabel('\beta^-')
name = join([dir,demo,'beta-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(5)
boxplot(lambda,table.age)
ylabel('\lambda')
name = join([dir,demo,'lambda.png']);
exportgraphics(gcf,name,'Resolution',300)
close

%% Errors box plot (panel reg)
errors = z_final(:,10);
lb = [0,0,0,0,0.01];
ub = [1,1,1,1,1];
cpt = z_final(:,5:9);
lb_error_norms = vecnorm(cpt - lb,2,2).^2;
ub_error_norms = vecnorm(cpt - ub,2,2).^2;
regs = z_final(:,4);
adjusted_errors = errors - regs.*lb_error_norms - regs.*ub_error_norms;

%%
% labels = {'$$\nu = 0$$','$$\nu = 1$$','$$\nu = 0$$','Case 4'}
boxplot(adjusted_errors,regs,'Notch','on');
xlabel('Regularization parameter $$\nu$$','Interpreter','latex');
ylabel('Squared norm of adjusted errors','Interpreter','latex');

%% Errors box plot (lottery)
errors = zeros(664,4);
errors(:,1) = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Code\AAC_cpt\Survey analysis\CPT\error_lottery_reg0_MODweight_cdf_norm-maxOutcome.mat').error;
errors(:,2) = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Code\AAC_cpt\Survey analysis\CPT\error_lottery_reg0_MODweight_cdf_norm-trueCE.mat').error;
errors(:,3) = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Code\AAC_cpt\Survey analysis\CPT\error_lottery_reg0_OGweight_cdf_norm-maxOutcome.mat').error;
errors(:,4) = load('C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\Code\AAC_cpt\Survey analysis\CPT\error_lottery_reg0_OGweight_cdf_norm-trueCE.mat').error;
boxplot(errors,'Notch','on','Labels',{'Case 1','Case 2','Case 3','Case 4'},'DataLim',[-Inf,5],'ExtremeMode','compress')
ylabel('Squared norm of normalized errors','Interpreter','latex');

