%% Loading data
cd 'C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC'
T = load('Data\Panel\Full_launch_6-10_n955');
table_orig = T.Fulllaunch610n955;
cd 'Code\AAC_cpt\Survey analysis'
table = convertvars(table_orig,[7:8,10:12,48:139],'double');

cpt = load('Code\AAC_cpt\Survey Analysis\Results\cpt_lottery_improved');
error = load('Code\AAC_cpt\Survey Analysis\Results\error_cpt_lottery_improved');
load('Code\AAC_cpt\Survey Analysis\valid_indices.mat')
[I,I_id] = findgroups(table.income);

%%
num_valid = length(valid_indices);
respondents = 1:1:num_valid;

%% Small scale with only 20 respondents
num = len(valid_indices);
respondents = 1:1:num;
cpt = cpt(1:num,:);
error = error(1:num,:);

alpha_plus = cpt(:,1);
alpha_minus = cpt(:,2);
beta_plus = cpt(:,3);
beta_minus = cpt(:,4);
lambda = cpt(:,5);

% Scatter plots

figure(1)
scatter(respondents,alpha_plus)
ylabel('\alpha^+')
xlabel('Respondent')

figure(2)
scatter(respondents,alpha_minus)
ylabel('\alpha^-')
xlabel('Respondent')

figure(3)
scatter(respondents,beta_plus)
ylabel('\beta^+')
xlabel('Respondent')

figure(4)
scatter(respondents,beta_minus)
ylabel('\beta^-')
xlabel('Respondent')

figure(5)
scatter(respondents,lambda)
ylabel('\lambda')
xlabel('Respondent')

figure(6)
scatter(respondents,error)
ylabel('Squared norm of residual')
xlabel('Respondent')

%% Histograms
figure(1)
histogram(alpha_plus);
xlabel('\alpha^+')
ylabel('Number of respondents');

figure(2)
histogram(alpha_minus);
xlabel('\alpha^-')
ylabel('Number of respondents');

figure(3)
histogram(beta_plus);
xlabel('\beta^+')
ylabel('Number of respondents');

figure(4)
histogram(beta_minus);
xlabel('\beta^-')
ylabel('Number of respondents');

figure(5)
histogram(lambda);
xlabel('\lambda')
ylabel('Number of respondents');

%% Population averages
% Mean
alpha_plus_mean = mean(alpha_plus)
alpha_minus_mean = mean(alpha_minus)
beta_plus_mean = mean(beta_plus)
beta_minus_mean = mean(beta_minus)
lambda_mean = mean(lambda)

%% Median
clc;
alpha_plus_median = median(alpha_plus)
alpha_minus_median = median(alpha_minus)
beta_plus_median = median(beta_plus)
beta_minus_median = median(beta_minus)
lambda_median = median(lambda)

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

%% Standard errors
clc;
alpha_plus_sd = std(alpha_plus)/sqrt(955)
alpha_minus_sd = std(alpha_minus)/sqrt(955)
beta_plus_sd = std(beta_plus)/sqrt(955)
beta_minus_sd = std(beta_minus)/sqrt(955)
lambda_sd = std(lambda)/sqrt(955)

%% Box plots 
% Gender
figure(1)
boxplot(alpha_plus,table.gender)
ylabel('\alpha^+')
xlabel('Gender');

figure(2)
boxplot(alpha_minus,table.gender)
ylabel('\alpha^-')
xlabel('Gender');

figure(3)
boxplot(beta_plus,table.gender)
ylabel('\beta^+')
xlabel('Gender');

figure(4)
boxplot(beta_minus,table.gender)
ylabel('\beta^-')
xlabel('Gender');

figure(5)
boxplot(lambda,table.gender)
ylabel('\lambda')
xlabel('Gender');

%% Income
Income_order = string([I_id(11); I_id(1); I_id(3:10); I_id(2); I_id(12)]);
figure(1)
boxplot(alpha_plus,table.income,'PlotStyle','compact','GroupOrder',Income_order)
ylabel('\alpha^+')
xlabel('Income');

figure(2)
boxplot(alpha_minus,table.income,'PlotStyle','compact','GroupOrder',Income_order)
ylabel('\alpha^-')
xlabel('Income');

figure(3)
boxplot(beta_plus,table.income,'PlotStyle','compact','GroupOrder',Income_order)
ylabel('\beta^+')
xlabel('Income');

figure(4)
boxplot(beta_minus,table.income,'PlotStyle','compact','GroupOrder',Income_order)
ylabel('\beta^-')
xlabel('Income');

figure(5)
boxplot(lambda,table.income,'PlotStyle','compact','GroupOrder',Income_order)
ylabel('\lambda')
xlabel('Income');

%% Age
figure(1)
boxplot(alpha_plus,table.age)
ylabel('\alpha^+')
xlabel('Age');

figure(2)
boxplot(alpha_minus,table.age)
ylabel('\alpha^-')
xlabel('Age');

figure(3)
boxplot(beta_plus,table.age)
ylabel('\beta^+')
xlabel('Age');

figure(4)
boxplot(beta_minus,table.age)
ylabel('\beta^-')
xlabel('Age');

figure(5)
boxplot(lambda,table.age)
ylabel('\lambda')
xlabel('Age');
