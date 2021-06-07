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

z = array2table(z_final(:,3:10),'VariableNames',...
    {'Respondent','nu','alpha_plus','alpha_minus','beta_plus','beta_minus','lambda','error'});

reg = 1.5;
rows = (round(z.nu,1) == reg);
cpt = table2array(z(rows,3:7));
error = table2array(z(rows,8));

alpha_plus = cpt(:,1);
alpha_minus = cpt(:,2);
beta_plus = cpt(:,3);
beta_minus = cpt(:,4);
lambda = cpt(:,5).*100;

%% Histograms
figure(1)
histogram(alpha_plus)
xlabel('$\alpha^+$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')
set(gca,'FontSize',14) 

figure(2)
histogram(alpha_minus)
xlabel('$\alpha^-$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')
set(gca,'FontSize',14) 

figure(3)
histogram(beta_plus)
xlabel('$\beta^+$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')
set(gca,'FontSize',14) 

figure(4)
histogram(beta_minus)
xlabel('$\beta^-$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')
set(gca,'FontSize',14) 

figure(5)
histogram(lambda)
xlabel('$\lambda$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')
set(gca,'FontSize',14) 

%%
non_outliers = lambda <= 20;
new_lambda = lambda(non_outliers);
histogram(new_lambda);
xlabel('$\lambda$','Interpreter','latex')
ylabel('Number of respondents','Interpreter','latex')
set(gca,'FontSize',14) 

%%
edges = [0:5:90 100];
h = histogram(lambda, edges);
h.Normalization = 'countdensity';