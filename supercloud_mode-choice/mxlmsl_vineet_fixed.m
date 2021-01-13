 % Matlab code to estimate a mixed logit model with maximum simulated likelihood
% For use on MIT Supercloud
% Need to store parameter estimates, standard errors and log-likelihood history 
% Loop over combinations of hyperparams (no. of draws, draw type, distribution type, starting values)
% Choose hyperparam that maximizes likelihood ratio index

clc;
clear;

% Declare GLOBAL variables
% GLOBAL variables are all in caps
% DO NOT CHANGE ANY OF THESE 'global' STATEMENTS
global NP NCS NROWS
global IDV NV NAMES B W
global IDF NF NAMESF F
global DRAWTYPE NDRAWS SEED1 SAVEDR PUTDR
global WANTWGT IDWGT WGT
global NALTMAX NCSMAX
global X XF S DR
global XMAT
global NMEM NTAKES NPARAM
global MDR

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

% delete myrun.txt
% diary myrun.txt

% TITLE
% Put a title for the run in the quotes below, to be printed at the top of the output file.

% DATA

% Number of people (decision-makers) in dataset 
NP=955;        

% Number of choice situations in dataset. This is the number faced by all the people combined.
NCS=NP*11;     

% Total number of alternatives faced by all people in all choice situations combined.
% This is the number of rows of data in XMAT below.
NROWS=NCS*3;

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

XMAT=load('Xmat.txt'); %The variables are described below

XMAT(:,4:9)= -XMAT(:,4:9); %To make price and travel times negative so coefficients can be positive.

% To help you keep up with the variables, list the variables in XMAT here.
% This dataset is for people's choice among vehicles in stated-preference
% experiments. Each person faced 11 experiments 
% Each experiment contained 3 alternatives representing different travel
% options whose price and other attributes were described. The person stated which
% of the 3 modes he/she would select if facing this choice in the real world.
% The variables in XMAT are:
% 1. Person number (1-NP)            MUST BE THIS. DO NOT CHANGE.
% 2. Choice situation number (1-NCS) MUST BE THIS. DO NOT CHANGE.
% 3. Chosen alternative (1/0)        MUST BE THIS. DO NOT CHANGE.
% 4. Walking  time
% 5. Waiting time
% 6. Transit riding time
% 7. Exclusive rideshare riding time
% 8. Pooled rideshare riding time
% 9. Price 
% 10. Public transit (1/0)
% 11. Exclusive rideshare (1/0)
% 12. Pooled rideshare (1/0)

% MODEL SPECIFICATION

% RANDOM COEFFICIENTS
% List the variables in XMAT that enter the model with random coefficients and
% give the distribution for the coefficient of each variable.
% IDV contains one row for each random coefficient and two columns.
% The *first* column gives the number of a variable in XMAT that has a random coefficient, 
% and the *second* column specifies the distribution of the coefficient for that variable.

% The distributions can be 
% 1. normal: N(b,w^2) where mean b and standard deviation w are estimated.
% 2. lognormal: coefficient is exp(beta) where beta~N(b,w^2) with b and w estimated
% 3. truncated normal, with the share below zero massed at zero: max(0,beta) where 
%                      beta~N(b,w^2) with b and w estimated.
% 4. S_B: exp(beta)/(1+exp(beta))  where beta~N(b,w^2) with b and w estimated.
% 5. normal with zero mean (for error components): N(0,w^2) where w is estimated.
% 6. triangular: b+w*t where t is triangular between -1 and 1 and mean b and spread w are estimated.

% If no random coefficients, put IDV=[];
% Notes:
% The lognormal, truncated normal, and S_B distributions give positive
% coefficients only. If you want a variable to have only negative coefficients, 
% create the negative of the variable (in the specification of XMAT above).
% The S_B distribution gives coefficients between 0 and 1. If you want
% coefficients to be between 0 and k, then multiply the variable by k (in the specification 
% of XMAT above), since b*k*x for b~(0-1) is the same as b*x for b~(0-k).
% If no random coefficients, put IDV=[];

IDV=[];

NV=size(IDV,1); %Number of random coefficients. Do not change this line.

% Give a name to each of the explanatory variables in IDV. They can 
% have up to ten characters including spaces. Put the names in single quotes and separate 
% the quotes with semicolons. If IDV=[], then set NAMES=[];
NAMES=[];

% Starting values
% Specify the starting values for b and w for each random coeffient.
% B contains the first parameter, b, for each random coefficient.  
% It is a column vector with the same length as IDV. For distribution 5 (normal with zero mean),
% put 0 for the starting value for the mean. The code will keep it at 0.
% W contains the second parameter, w, for each random coefficient.
% It is a column vector with the same length as IDV.
% Put semicolons between the elements of B and W (so they will be column vectors).

% Set both of these to be zero? So that LL at 1st iteration is with all model params = 0
% Or do I need to calculate LL(0) using all fixed coefficients with starting value 0?
% Then use this to calculate likelihood ratio index 

B=[];
W=[];

% FIXED COEFFICIENTS
% List the variables in XMAT that enter with fixed coefficients.
% Put semicolons between the numbers.
% If no fixed coefficients, put IDF=[];

% Set public transit as the baseline alternative with ASC_transit = 0
IDF=[4;5;6;7;8;9;11;12]; 

NF=size(IDF,1); %Number of fixed coefficients. Do not change this line.

% Give a name to each of the variables in IDF.
NAMESF={'t_walk';'t_wait';'tr_ride';'exc_rideshare';'pool_rideshare';'price';'exclusive';'pooled'}; 

% Starting values.
% Specify the starting values for the fixed coefficients F.
% F must have the same length as IDF and have one column.
% Put semicolons between the elements (so F will be a column vector.)

F=[init; init; init; init; init; init; init; init];

% Type of draws to use in simulation
% 1=pseudo-random draws
% 2=standard Halton draws
% 3=shifted and shuffled Halton draws
% 4=modified Latin hypercube sampling, shifted and shuffled 
% 5=create your own draws or load draws from file
DRAWTYPE=draw_type;

% Number of draws from to use per person in simulation.
NDRAWS=ndraws;

% Set seed for the random number generator.
SEED1 = 14239; 

% If DRAWTYPE=5, then create or load the draws here.
% Create or load a data array, called DR, with dimensions NV x NP x NDRAWS.
% where element DR(i,j,k) is the k-th draw of random coefficient i for person j 
% Put your statements between "if DRAWTYPE==5" and "end".

% If you created DR in a previous matlab session and saved it
% with "save mydraws DR" then put "load('mydraws.mat')" here. The structure
% mydraws will contain the array DR.
% Note: If you want to use the draws that were saved to PUTDR below in a 
% previous run, see the ReadMe.txt file for instructions.

if DRAWTYPE==5
   load('mydraws.mat');
end

% Memory use
% Give the number of draws that you want held in memory at one time.
% This number must be evenly divisible into the number of draws.
% That is NDRAWS./NMEM must be a positive integer.
% To hold all draws in memory at once, set NMEM=NDRAWS.
% A larger value of NMEM requires fewer reads from disc but 
% uses more memory which can slow-down the calculations and increases 
% the chance of running out of memory.
% If DRAWTYPE=5, then you must set NMEM=NDRAWS
NMEM=NDRAWS;

% If all the draws are NOT held in memory at one time (that is, if NMEM<NDRAWS), 
% then give the filename (including full path if not in the working directory)
% that you want the draws to be temporarily saved to while the code is running.
% If all draws are held in memory at one time (that is, if NMEM=NDRAWS),
% then this file will not be created. So, if NMEM=NDRAWS, you can set PUTDR=''; 
% or give a file name, whichever you find more convenient, since the name won't be used.
PUTDR='draws';

% WEIGHTS. 
% Do you want to apply weights to the people? 
% Set WANTWGT=1 if you want to apply weights; otherwise set WANTWGT=0;
WANTWGT=0;

% If WANTWGT=1, identify the variable in XMAT that contains the weights.
% This variable can vary over people but must be the same for all rows of
% data for each person. Weights cannot vary over choice situations for
% each person or over alternatives for each choice situation -- only over people.
% The code normalizes the weights such that the sum 
% of weights over people is to equal NP (to assure that standard errors 
% are correctly calculated.) If WANTWGT=0, set IDWGT=[];
IDWGT=[];

% OPTIMIZATION 
% Maximum number of iterations for the optimization routine.
% The code will abort after ITERMAX iterations, even if convergence has
% not been achieved. The default is 400, which is used when MAXITERS=[];
MAXITERS=1000;

% Convergence criterion based on the maximum change in parameters that is considered
% to represent convergence. If all the parameters change by less than PARAMTOL 
% from one iteration to the next, then the code considers convergence to have been
% achieved. The default is 0.000001, which is used when PARAMTOL=[];
PARAMTOL=[];

% Convergence criterion based on change in the log-likelihood that is
% considered to represent convergence. If the log-likelihood value changes
% less than LLTOL from one iteration to the next, then the optimization routine
% considers convergence to have been achieved. The default is 0.000001,
% which is used when LLTOL=[];
LLTOL=[];

%Do not change the next line. It runs the model.
doit
% These last lines delete the file of draws that is created when NMEM<NDRAWS
% since it is no longer needed. If you want to save it, then put % in front
% of these lines.
if NMEM<NDRAWS
    clear global MDR
    delete(PUTDR)
end
