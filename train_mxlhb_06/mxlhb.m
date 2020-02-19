% Matlab code to estimate a mixed logit model with hierachical Bayes procedures
% Written by Kenneth Train, first version July 13, 2006, 
%   latest edits July 27, 2006
%
% MODEL NOTATION
% Utility of person n from alternative j is: U_jn=F'Z_jn+C_n'X_jn+e_jn
% where F is a vector of fixed coefficients (the same for all people) of variables Z_jn. 
% and C_n is a vector of random coefficients (different for each person) of variables X_jn. 
% The random coefficients for each person are transformations of randomly distributed 
% latent terms. That is, C_n=T(B_n) for transformation T and B_n~N(A,D). 
% The transformations available in this code are:
%     normally distributed coefficient: C_n=B_n, 
%     lognormal: C_n=exp(B_n), 
%     truncated normal: C_n=max(0,B_n).
% The latent terms are distributed normally with mean vector A and covariance matrix D.
% Correlation among the latent terms induces correlation among the coefficients.
% The procedure takes draws of A, D, and B_n for all n sequentially using MCMC methods. 
% Draws of C_n are created from the draws of B_n. That is, a draw of B_n implies a draw of C_n.
% Matrix B={B_1,...,B_N} is a NVxNP matrix of latent terms for all people. 
% Matrix C is the corresponding matrix of coefficients for all people.

% DATA and MODEL SPECIFICATION

clear all

% Declare GLOBAL variables
% GLOBAL variables are all in caps
% DO NOT CHANGE ANY OF THESE 'global' STATEMENTS
global NP NCS NROWS
global IDV NV NAMES FULLCV
global A B C D  
global RHO NCREP NEREP NSKIP NRFOR NRLL INFOSKIP
global SEED1 SEED2
global KEEPMN PUTA PUTD WANTINDC PUTC 
global NALTMAX NCSMAX
global X S A B C D SQRTNP 
global XF F NF IDF NAMESF RHOF
global XMAT


% OUTPUT FILE
% Put the name you want for your output file (including full path if not the current 
% working directory) after words "delete" and "diary".
% The 'diary off' and 'delete filename' commands close and delete the previous version 
% of the file created during your current matlab session or in any previous sessions. 
% If you want to append the new output to the old output, then 
% put % in front of the 'diary off' and 'delete filename' commands (or erase them).

diary off
delete myrun.out
diary myrun.out

% TITLE
% Put a title for the run in the quotes below, to be printed at the top of the output file.
disp 'Mixed logit of vehicle choice.'


% DATA

% Number of people (decision-makers) in dataset 
NP=100;        

% Number of choice situations in dataset. This is the number faced by all the people combined.
NCS=1484;     

% Total number of alternatives faced by all people in all choice situations combined.
% This is the number of rows of data in XMAT below.
NROWS=1484.*3;

% Load and/or create XMAT, a matrix that contains the data.
%
% XMAT must contain one row of data for each alternative in each choice situation for each person.
% The rows are grouped by person, and by choice situations faced by each person.
% The number of rows in XMAT must be NROWS, specified above.
% The columns in XMAT are variable that describe the alternative.
% 
% The *first* column of XMAT identifies the person who faced this alternative. 
% The people must be numbered sequentially from 1 to NP, in ascending order.
% All alternatives for a given person must be grouped together.
% The *second* column of XMAT identifies the choice situation. The choice
% situations must be numbered sequentially from 1 to NCS.
% All alternatives for a given choice situation must be grouped together.
% The *third* column of XMAT identifies the chosen alternatives (1 for
% chosen, 0 for not). One and only one alternative must be chosen for each
% choice situation.
% The remaining columns of XMAT can be any variables.

XMAT=load('data.txt');  %The variables are described below

XMAT(:,4:5)= -XMAT(:,4:5);    %To make price and opcost negative so coef can be positive.
XMAT(:,4)=XMAT(:,4)./10000;   %To scale price to be in tens of thousands of dollars.

% To help you keep up with the variables, list the variables in XMAT here.
% Start each line with % so that matlab sees that it is a comment rather than a command.
% NOTES for XMAT for sample run:
% This dataset is for people's choice among vehicles in stated-preference
% experiments. Each person faced up to 15 experiments (some faced fewer
% than 15 because they did not complete all the experiments.) Each
% experiment contained 3 alternatives representing three different vehicles
% whose price and other attributes were described. The person stated which
% of the three vehicle he/she would buy if facing this choice in the real world.
% The variables in XMAT are:
% 1. Person number (1-NP)            MUST BE THIS. DO NOT CHANGE.
% 2. Choice situation number (1-NCS) MUST BE THIS. DO NOT CHANGE.
% 3. Chosen alternative (1/0)        MUST BE THIS. DO NOT CHANGE.
% 4. Negative of Price in tens of thousands of dollars 
% 5. Negative of Operating cost in dollars per month
% 6. Range in hundreds of miles (0 if not electric)
% 7. Electric (1/0)
% 8. Gas (1/0)
% 9. Hybrid (1/0)
% 10. High performance (1/0)
% 11. Medium or high performance (1/0)

% MODEL SPECIFICATION

% RANDOM COEFFICIENTS
% List the variables in XMAT that enter the model with random coefficients and
% give the distribution for the coefficient of each variable.
% IDV contains one row for each random coefficient and two columns.
% The *first* column gives the number of a variable in XMAT that has a random coefficient, 
% and the *second* column specifies the distribution of the coefficient for that variable.
% The distributions can be: 1. normal, 2. lognormal, 3. truncated normal (with the share 
% below zero massed at zero), 4. S_B, 5. normal with zero mean(for error components.)
% If no random coefficients, put IDV=[];
% Notes:
% The lognormal, truncated normal, and S_B distributions give positive
% coefficients only. If you want a variable to have only negative coefficients, 
% create the negative of the variable (in the specification of XMAT above).
% The S_B distribution gives coefficients between 0 and 1. If you want
% coefficients to be between 0 and k, then multiply the variable by k (in the specification 
% of XMAT above), since b*k*x for b~(0-1) is the same as b*x for b~(0-k).
 
IDV=[ 4 2;   ...
      5 2;   ...
      6 2;   ...
      7 1;   ...
      9 1];

NV=size(IDV,1); %Number of random coefficients. Do not change this line.

% Give a name to each of the explanatory variables in IDV. They can 
% have up to ten characters including spaces. Put the names in single quotes and separate 
% the quotes with semicolons.
NAMES={'price';'opcost';'range';'ev';'hybrid'}; 

% Correlation
% Set FULLCV=1 to estimate covariance among all random coefficients, FULLCV=0 for no covariances
FULLCV=1;

% Starting values
% Specify the starting values for A,B, and D.
% A contains the means of the distribution of the underlying normal for each random coefficient.  
% It is a column vector with the same length as IDV. For distribution 5 (normal with zero mean),
% put 0 for the starting value for the mean. The code will keep it at 0.
% B contains the random coefficients for each person.
% It is a matrix with IDV rows and NP columns.
% D contains the covariance matrix of the distribution of underlying normals for the 
% random coefficients. It is a symmetric matrix with IDV rows and columns.
% If FULLCV=0, then D is diagonal.

A=[ .1; .1;  .2 ; -1.44; .412];
B=repmat(A,1,NP);    %Start each person at the starting means
D=NV.*eye(NV);       %Starting variances equal to number of random coefficients

% Set the initial proportionality fraction for the jumping distribution for
% the random coefficients for each person.
% This fraction is adjusted by the program in each iteration to attain
% an acceptance rate of about .3 in the M-H algorithm for the B's 
RHO = 0.1;

% FIXED COEFFICIENTS
% List the variables in XMAT that enter with fixed coefficients.
% Put semicolons between the numbers.
% If no fixed coefficients, put IDF=[];

IDF=[10;11];

NF=size(IDF,1); %Number of fixed coefficients. Do not change this line.

% Give a name to each of the variables in IDF.
NAMESF={'hiperf'; 'medhiperf'};


% Starting values.
% Specify the starting values for the fixed coefficients F.
% F must have the same length as IDF and have one column.
F=[0;0];



% Set the proportionality fraction for the jumping distribution for
% the fixed coefficients.
% This fraction is not adjusted by the program. Ideally, RHOF should be set such
% the acceptance rate is about .3 in the M-H algorithm for the draws of F 
RHOF = 0.01;

% BURN-IN, RETAINING DRAWS, PRINTING ITERATION INFORMATION

% Number of iterations to make prior to retaining draws (i.e., length of burn-in)
NCREP = 10000;

% Number of draws to retain after burn-in 
NEREP = 1000;

% Number of iterations to make between retained draws
% NOTE: The number of iterations after burn-in is NEREP*NSKIP, of which NEREP are retained
% The total number of iterations is NCREP+(NEREP*NSKIP)
NSKIP = 10;

% Number of draws from N(A-hat,D-hat) to use in simulating the estimated
% dist of random coefficients.
NRFOR=2000;

% Number of draws per person to use in calculating the simulated log-likelihood value
% at A-hat, D-hat, and F-hat.
NRLL=2000;

% How frequently to print information about the iteration process 
% eg 100 == print out info every 100-th iteration 
INFOSKIP = 500;

% Set seeds for the random number generators.
% Seed for random draws in iteration process. Must be a positive integer.
SEED1 = 14239;   
% Seed for random draws fron N(A-hat,D-hat) in simulated distribution of
% coefficients and log-likelihood value. Must be a positive integer.
SEED2 = 123;   

% Do you want to save the draws of A, D, and F to a file? 
% If so, set KEEPMN=1. If no, set KEEPMN=0.
KEEPMN=1;

% If KEEPMN=1, specify filenames to which to save the draws of A, D, and F.
% Put the file name in single quotes.
% If you want the file saved to a directory other than the working
% directory, then put the entire path and filename (eg, 'c:\myfolder\drawsA').
% The program will save the NEREP draws of A to FOUTPUTA as a matlab matrix
% with dimension NVxNEREP
% The program will save the NEREP draws of D to FOUTPUTD as a  matlab array 
% of dimension NVxNVxNEREP.
% The program will save the NEREP draws of F to FOUTPUTF as a matlab matrix
% with dimension NFxNEREP.
% If the model does not contain random coefficients, then draws of A and D
% will not be saved (since they do not exist). Similarly, if the model does
% not contain fixed coefficients, then draws of F will not be saved.

PUTA = 'drawsA';
PUTD = 'drawsD';
PUTF = 'drawsF';

% Do you want to calculate and save the means of the NEREP draws of the
% individual-level random coefficients? If so, set WANTINDC=1. If no, set WANTINDC=0.
WANTINDC=1;

% If WANTINDC=1, specify filename (and full path if not the working directory) 
% to which to save the mean of the draws
% of the individual-level random coefficients. The output in FOUTPUTC
% will be a matlab matrix of dimension NPxNIV, containing one row for each
% person and one column for each coefficient. 
PUTC = 'meansC';

% Do not change the next line. It calls the script that performs the
% iterations and prints results.
doit
