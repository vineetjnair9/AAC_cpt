
param=StartB;

disp('Start estimation');
disp('The negative of the log-likelihood is minimized,');
disp('which is the same as maximizing the log-likelihood.');

options=optimset('LargeScale','off','Display','iter','GradObj','on',...
    'MaxFunEvals',10000,'MaxIter',2000,'TolX',10^(-6),'TolFun',10^(-6),'DerivativeCheck','off');

[paramhat,fval,exitflag]=fminunc(@flexll,param,options);

disp('Calculating summary statistics for random utility parameters.');
load('test_save.mat');
[MeanEst,StdEst,CovMatEst,FreqEst,MidEst]=stats(paramhat,NBins);
disp('Estimated coefficients of Z variables and non-random are held in paramhat.');

disp('Means, Standard Deviations, and Covariance Matrix of Utility Parameters');
disp('are held as MeanEst, StdEst, and CovMatEst.');
disp('Share of density in each bin and midpoint for each bin are held');
disp('in FreqEst and MidEst; each is size NV x NBins');
disp('Means and StdDevs of Random Utility Paramaters');
disp(' ');
disp('                Mean');
disp('              -------------------');
disp('              Estimate      True');
for r=1:length(NAMES)
    fprintf('%-10s %10.4f %10.4f\n', NAMES{r,1}, MeanEst(r,1), true_par_mean(r,1));
end


disp(' ');
disp('Var-covar Matrix for Random Coefficients');
disp(CovMatEst);

disp(' ');
disp('True Var-covar Matrix for Random Coefficients');
disp(true_covar);

disp(' ');
disp('To get histogram of random coefficient 1');
disp('use bar(MidEst(1,:),FreqEst(1,:))');


