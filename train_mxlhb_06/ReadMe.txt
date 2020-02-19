Mixed Logit Estimation by Hierarchical Bayes, written in Matlab by Kenneth Train. Current version dated September 24, 2006.

The estimation procedure is described in:
Kenneth Train, Discrete Choice Methods with Simulation, Chapter 12 (New York: Cambridge University Press, 2003.) 
and
Kenneth Train and Garrett Sonnier, "Mixed Logit with Bounded Distributions of Correlated Partworths,"
Chapter 7 of Applications of Simulation Methods in Environmental and Resource Economics, Riccardo
Scarpa and Anna Alberini, editors (Dordrecht, The Netherlands: Springer, 2005.)

Both publications are available on-line at http://elsa.berkeley.edu/~train.

Notation within codes: F is vector of fixed coefficients. C is a matrix of random coefficients (a vector for each person). The random coefficients are transformations of latent random terms B that are normally distributed with mean A and covariance D.

Put all the files into one directory and use that directory when running Matlab (or change paths as needed.) To check that everything is working correctly, run mxlhb.m as is. The output should be the same as in myrunKT.out (unless you are using a very early version of matlab that has a different random number generators, in which case the results should be similar but not exactly the same.)

FILES:

mxlhb.m is the code that the user runs. The user specifies the model within this code. 

doit.m is a script (not a function) that is called at the end of mxlhb.m. It checks the data,
transforms the data into a more useful form, performs the estimation and prints results. It calls all the other functions either directly or indirectly.

check.m is a function that checks the input data and specifications. It provides error messages and
terminates the run if anything is found to be incorrect.

hbmcmc.m is a function that performs the iterative mcmc procedure and prints the iteration log,
calling the functions below.

trans.m is a function that transforms the latent normal terms B into coefficients C, using the user's specifications in mxlhb.m. 

logit.m is a function that calculates the logit probaility for each person, using fixed 
coefficients F and random coefficients C as input arguments.

nextB.m is a function that calculates the a draw of B, the latent normal terms for each person, 
using an MH algorithm on each person separately. A draw of B implies a draw of C through trans.m.

nextA.m is a function that calculates a draw of A, the mean of the latent normal terms.

nextD.m is a function that calculates a draw of D, the covariance matrix of the latent normal terms,
given that D is full (ie allows correlation among all latent terms.)

nextDdiag.m is a function that calculates a draw of D, the covariance matrix of the latent normal
terms, given that D is specified by the user to be diagonal (no correlations.)

nextF.m is a function that calculates a draw of F, the fixed coefficients that enter utility, 
using an MH algorithm on entire sample.

data.txt is an ascii file of data on vehicle choice. The data and its format are described
within mxlhb.m

myrunKT.out is the output file of running mxlhb.m with no modifications (with KT added to the
filename so it is not deleted when the code is run again.) 

FILE HIERARCHY

mxlhb calls doit. mxlhb also reads in data.txt and outputs the diary of the run and, if requested,
the retained draws of A, D, and F and the person-specific means of the draws of C.

doit calls check to check the data and then hbmcmc to perform the iterations. It also calls 
trans and logit to calculate the log-likelihood value at mean of retained draws of A,D, and F.

hbmcmc calls nextB, nextA, nextD, nextDdiag, nextF iteratively to perform the mcmc procedure. It
calls trans and logit to start the process.

nextB and nextF call trans and logit. 
