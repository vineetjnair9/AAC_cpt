% Matlab code to estimate a flexible mixed logit model in PREFERENCE space(with random and non-random coefficient).
% Written by Kenneth Train, first version Oct 21, 2015, Modified by Prateek Bansal, on Feb 23, 2019.

clear all

% Do not change the next line. It sets the global variables.
global NP NCS NROWS NZ  NDRAWS NV IDV ZTYPE CrossCorr COEF NNR IDNR NAMESNR XMAT NAMESV NAMES NBins BETAS Z NReps

% Load and/or create XMAT, a matrix that contains the data.
XMAT = csvread('test.csv');

% DATA
% Number of people (decision-makers) in dataset
NP=500;

% Number of choice situations in dataset. This is the number faced by all the people combined.
NCS=NP*8;

% Total number of alternatives faced by all people in all choice situations combined.
NROWS=NCS*3;


% The *first* column of XMAT identifies the person who faced this alternative.
% The people must be numbered sequentially from 1 to NP, in ascending order.
% All alternatives for a given person must be grouped together.
% The *second* column of XMAT identifies the choice situation. The choice
% situations must be numbered sequentially from 1 to NCS.
% All alternatives for a given choice situation must be grouped together.
% The *third* column of XMAT identifies the chosen alternatives (1 for
% chosen, 0 for not). One and only one alternative must be chosen for each
% choice situation.

% MODEL SPECIFICATION

% RANDOM COEFFICIENTS
% Put semicolons between the numbers (so that IDV is a column vector).
IDV = [6;7];% Column ID of random coefficient
IDNR = [4;5];% Column ID of non-random coefficients


NV=size(IDV,1);
NNR = size(IDNR,1);

NAMESV={'random1';'random2'};
NAMESNR = {'fixed1';'fixed2'};
NAMES = {'random1';'random2';'fixed1';'fixed2'};


% Give the range for the coeff for each variable with RANDOM coefficient.
% The first column gives the lower limit and the second column gives the upper limit.
% Limits are inclusive, in that value of the limit is included.
% Put colons between range for each variable.
COEF =[-4    2;
    -2  4];  

%Give the number of points in each dimension to define the grid for the
%parameter space. The grid includes the two endpoints given in COEF.
%So if you want a grid with eg 1000 intervals between these endpoints, then setNGridPts=1001.
%The total number of points in the grid is NGridPts^NV
NGridPts=1000;

%Number of random draws from the parameter space to use in simulation for each person
NDRAWS=2000;

%Specify the seed to use in the random number generator for simulation
ThisSeed=1234;

%No not change the next line
rng(ThisSeed);

% Specify the type of variables that you want to describe the marginal densities of
% the random parameters. The options are:
% ZTYPE=1 for polynomials.
%      =2 for step functions
%      =3 for splines
ZTYPE=1;

%If ZTYPE=1, specify the order of the polynomial.
%The number of Z variables that will be created is PolyOrder*NV where
%where NV is the number of random utility parameters.
PolyOrder=2;

%If ZTYPE=2, specify the number of levels for the step functions.
%The height of each level is estimated, with the height of the final step
%normalized to zero. The number of Z variables that will be
%created is (NLevels-1)*NV
NLevels=3;

%If ZTYPE=3, specify the number of knots for each spline.
%The spline starts at the low end of the range given above,
%changes slopes at each knot, and stops at the high end of the range.
%The number of parameters for a spline with K knots is K+1, i.e.,
%the height at the low end and each knot, with the height
%at the high end normalized to zero. The number of Z variables that
%will be created is (NKnots+1)*NV
NKnots=2;

%Do you want to include cross-products for correlations among utility coefficients?
%Set CrossCorr=1 for yes, =0 for no.
%If ZTYPE=1, then the number of extra Z variables that are created with
%CrossCorr=1 is (NV)*(NV-1)/2. (There are NV coefficients, and so
%(NV)*(NV-1)/2 pairs of coefficients.)
%   If ZTYPE=2 or 3, then the number of extra Z variables is
% 2*(NV)+[(NV)*(NV-1)/2], i.e, a second order polynomial in each random coefficient by
% itself (for 2*(NV) extra Z variables) plus cross-product terms for each pair
% of random coefficients.

CrossCorr=1;

%Do not change the following lines. They calculate the number of Z variables
if ZTYPE==1
    NZ=PolyOrder*NV;
    if CrossCorr==1
        NZ=NZ+(NV)*(NV-1)/2;
    end
elseif ZTYPE==2
    NZ=(NLevels-1)*NV;
    if CrossCorr==1
        NZ=NZ+2*(NV)+(NV)*(NV-1)/2;
    end
elseif ZTYPE==3
    NZ=(NKnots+1)*NV;
    if CrossCorr==1
        NZ=NZ+2*(NV)+(NV)*(NV-1)/2;
    end
end

%Set the starting values for the coefficients of the Z variables and
%non-random coefficient
StartB=zeros(NZ,1);  % starting value of z variables
StartNR = zeros(NNR,1); % starting value of non-random coefficient
StartB = [StartB; StartNR];
%The code will estimate the coefficients of the Z variables and create
% a histogram for each random utility parameter based on the estimated
% coefficients of the Z variables.
%Specify the number of bins you want in the histogram.
%More bins gives more detailed shapes but results in more simulation noise in each bin.
NBins=20;

%Do you want to bootstrap the standard errors for estimated coefficients?
%Set WantBoot=1 for yes, =0 for no.
%Bootstrapping takes much longer than original estimation, and so
%you might want to bootstrap only after you are pretty sure of your model
%If WantBoot=1, specify the number of resamples:
WantBoot=0;
NReps=5;

DrawsBeta;
CreateZ;
DoEstimation;

if WantBoot==1
    Boot
end