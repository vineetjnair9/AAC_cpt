clc;
clear;
rng(10)
%% data generation for (Mixed Mixed Logit with fixed latent class coefficients)
N = 500; %the number of individuals
T = 8;  % the number of time periods observed for each individual
nrx = 2; % number of x variables with random coefficients
nfx =2; % number of observed variables with fixed coefficients
J =3; % the number of choices

FC  = [-0.8;0.8];
RCmean = [-1;1]; 
RCvar = [1,.7; .7, 1]; 

NROWS = N*T*J;
NCS = N*T;

NCOLS = 4 + nfx  + nrx;
data = zeros(NROWS,NCOLS);

%%Personid
col1 = repmat(1:N,J*T,1);
data(:,1) = reshape(col1,NROWS,1);
%choice situation id
col2 = repmat(1:NCS,J,1);
data(:,2) = reshape(col2,NROWS,1);
%col3 = chosen alternative
data(:,4:(NCOLS-1)) = normrnd(0,1,[NROWS (nfx+nrx)]); %covariates (first 2 fixed, then 2 random)

col8 =repmat(1:J,NCS,1);
data(:,NCOLS) = reshape(col8',NROWS,1);  %alternative numbering
%% generating datasets
temp_par = [repmat(FC,1,N)' mvnrnd(RCmean,RCvar,N) ];
parameters = kron(temp_par,ones(T*J,1));


error=evrnd(0,1,[NROWS,1]); %Different realizations of error term for each dataset
utility = sum(data(:,4:(NCOLS-1)).*parameters,2) + error;
V = sum(data(:,4:(NCOLS-1)).*parameters,2);
j=1;
keep_ind = zeros(NCS,1);

k=1;
counterror = 0;
while j <= NROWS
    [value,ind] = max(utility(j:(j+J-1)));
    [value1,ind1] = max(V(j:(j+J-1)));
    if ind1 ~= ind
        counterror = counterror+1;
    end
    
    data(j+ind-1,3)=1;
    keep_ind(k) = ind;
    j=j+J;
    k=k+1;
end

tabulate(keep_ind);
error_percent = (counterror/(N*T))*100;
csvwrite('test.csv',data)

true_covar = RCvar;
true_par_mean = [RCmean;FC];
save('test_save.mat','true_covar','true_par_mean','temp_par')


