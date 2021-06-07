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

cpt = cpt_lottery;
alpha_plus = cpt(:,1);
alpha_minus = cpt(:,2);
beta_plus = cpt(:,3);
beta_minus = cpt(:,4);
lambda = cpt(:,5);

%% Convert z_final array to table
z = array2table(z_final(:,3:10),'VariableNames',...
    {'Respondent','nu','alpha_plus','alpha_minus','beta_plus','beta_minus','lambda','error'});

clear T z_final table_orig

%% Swarm charts
% Need to change to nicer colormap!
dir = 'C:\Users\vinee\OneDrive - Massachusetts Institute of Technology\MIT\AAC\SM Thesis\Figures\CPT_demographic\';
demo = 'gender';

colormap winter
% Gender
figure(1)
swarmchart(table.gender,alpha_plus,20,grp2idx(table.gender),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\alpha^+')
name = join([dir,demo,'swarm_alpha+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(2)
swarmchart(table.gender,alpha_minus,20,grp2idx(table.gender),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\alpha^-')
name = join([dir,demo,'swarm_alpha-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(3)
swarmchart(table.gender,beta_plus,20,grp2idx(table.gender),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\beta^+')
name = join([dir,demo,'swarm_beta^+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(4)
swarmchart(table.gender,beta_minus,20,grp2idx(table.gender),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
name = join([dir,demo,'swarm_beta-.png']);
ylabel('\beta^-')
exportgraphics(gcf,name,'Resolution',300)
close

figure(5)
swarmchart(table.gender,lambda,20,grp2idx(table.gender),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
name = join([dir,demo,'swarm_lambda.png']);
ylabel('\lambda')
exportgraphics(gcf,name,'Resolution',300)
close

%% Income
demo = 'income';
% income_groups = categories(table.income);
income_groups = renamecats(table.income,{'$10-20K','$100-150K','$20-30K','$30-40K','$40-50K','$50-60K','$60-70K','$70-80K','$80-90K','$90-100K','< $10K','>$150K'});
income_groups = reordercats(income_groups,{'< $10K','$10-20K','$20-30K','$30-40K','$40-50K','$50-60K','$60-70K','$70-80K','$80-90K','$90-100K','$100-150K','>$150K'});
%%
figure(1)
swarmchart(income_groups,alpha_plus,20,grp2idx(income_groups),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\alpha^+')
name = join([dir,demo,'swarm_alpha+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(2)
swarmchart(income_groups,alpha_minus,20,grp2idx(income_groups),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\alpha^-')
name = join([dir,demo,'swarm_alpha-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(3)
swarmchart(income_groups,beta_plus,20,grp2idx(income_groups),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\beta^+')
name = join([dir,demo,'swarm_beta+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(4)
swarmchart(income_groups,beta_minus,20,grp2idx(income_groups),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\beta^-')
name = join([dir,demo,'swarm_beta^-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(5)
swarmchart(income_groups,lambda,20,grp2idx(income_groups),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\lambda')
name = join([dir,demo,'swarm_lambda.png']);
exportgraphics(gcf,name,'Resolution',300)
close

%% Age
demo = 'age';
figure(1)
swarmchart(table.age,alpha_plus,20,grp2idx(table.age),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\alpha^+')
name = join([dir,demo,'swarm_alpha+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(2)
swarmchart(table.age,alpha_minus,20,grp2idx(table.age),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\alpha^-')
name = join([dir,demo,'swarm_alpha-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(3)
swarmchart(table.age,beta_plus,20,grp2idx(table.age),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\beta^+')
name = join([dir,demo,'swarm_beta+.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(4)
swarmchart(table.age,beta_minus,20,grp2idx(table.age),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\beta^-')
name = join([dir,demo,'swarm_beta^-.png']);
exportgraphics(gcf,name,'Resolution',300)
close

figure(5)
swarmchart(table.age,lambda,20,grp2idx(table.age),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5)
ylabel('\lambda')
name = join([dir,demo,'swarm_lambda.png']);
exportgraphics(gcf,name,'Resolution',300)
close

%% Correlations
cpt = cpt_lottery;
corrplot(cpt,'varNames',{'a','b','c','d','e'})


