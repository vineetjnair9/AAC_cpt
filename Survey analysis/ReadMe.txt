Mixed Logit Estimation by Maximum Simulated Likelihood, written in Matlab by Kenneth Train. Current version dated September 24, 2006.

The model and estimation procedure are described in:
David Revelt and Kenneth Train, "Mixed Logit with Repeated Choices: Households' Choices of Appliance Efficiency Level," Review of Economics and Statistics, Vol. 80, No. 4, Nov. 1998, pp. 647-657.
and
Kenneth Train, "Recreation Demand Models with Taste Differences Over People," Land Economics, Vol. 74, No. 2, May 1998, pp. 230-239.
and 
Kenneth Train, Discrete Choice Methods with Simulation, Chapters 6 and 10 (New York: Cambridge University Press, 2003.) 

The data are a subset of the data described in:
Kenneth Train and Garrett Sonnier, "Mixed Logit with Bounded Distributions of Correlated Partworths," Ch. 7 (pp. 117-134) in Applications of Simulation Methods in Environmental and Resource Economics, Riccardo Scarpa and Anna Alberini, editors, (Dordrecht, The Netherlands: Springer, 2005).

All of these publications are available on-line at http://elsa.berkeley.edu/~train.

To run the code, you need to have Matlab and the Matlab Optimization Toolbox installed on your machine. Put all the files from mxlmsl.zip into a directory and use that directory as your working directory in Matlab (or change paths as needed.) To check that everything is working correctly, run mxlmsl.m. The output should be the same as in myrunKT.out (unless you are using a very early version of matlab that has a different random number generators, in which case the results should be similar but not exactly the same.)

FILES:

mxlmsl.m is the code that the user runs. The user specifies the model within this code. 

doit.m is a script (not a function) that is called at the end of mxlmsl.m. It checks the data,
transforms the data into a more useful form, performs the estimation and prints results. It calls all the other functions either directly or indirectly.

check.m is a function that checks the input data and specifications. It provides error messages and
terminates the run if anything is found to be incorrect.

loglik.m is a function that calculates the log-likehood function and its gradient. This funtion is input to Matlab's fminunc command (which is part of Matlab's Optimization Toolbox.) This function calls llgrad2.m.
 
llgrad2.m is a function that calculates for each person the probability of the chosen alternatives and the gradient of the log of this probability. 

der.m is a function that calculates the derivative of each random coefficient with respect to the parameters of the model.

makedraws.m is a function that creates the standardized (ie parameter-free) draws that will be used in the run, based on the specifications given by the user in mxlmsl.m. 

trans.m is a function that transforms the standardized draws into draws of coefficients. (EG, if coefficient c is normal with mean b and standard error w, then makedraws.m creates draws mu from a standard normal N(0,1), and trans.m creates the draws of coefficients as c=b+w*mu.)

data.txt is an ascii file of data on vehicle choice. The data and its format are described
within mxlhb.m

myrunKT.out is the output file of running maxlhb.m with no modifications (with KT added to the
filename so it is not deleted when the code is run again.) 

FILE HIERARCHY

mxlmsl calls doit. mxlmsl also reads in data.txt and outputs the diary of the run.

doit calls check to check the data, calls makedraws to make the standardized draws, and calls loglik as an input to matlab command fminunc to perform the estimation.

loglik calls llgrad2, which calls trans and der


NOTE: The file mxlmsl.m provides the user with an option to hold in memory only a portion of the draws that are used in estimation. This capability allows the users with large-scale models (many people and/or many random coefficients) to specify more draws than could be held in memory at one time without running out of memory. If this option is chosen, then makedraws.m saves the standardized draws to a file specified by the user. Then when calculating probabilities and gradients, llgrad2 "loops" over portions of this file, holding only part of it in memory each time. The file is deleted at the end of the run. If you want to access the file for some reason, then remove (or comment-out) the line "delete(PUTDR)" at the end of mxlmsl. The file is created with Matlab's fwrite, which means that it is not interpretable in a text editor nor can it be loaded into Matlab with a load command. If you want to see or use the contents, you can read the entire file into Matlab with these commands:

m=memmapfile('filename','Format','double');
DR=m.Data;
DR=reshape(DR,NDRAWS,NP,NV);
DR=permute(DR,[3,2,1]);

where filename is the name you specified for the file with PUTDR='filename' in the run that created it;
NDRAWS is the number of draws that you specified in the run that created the file.
NP is the number of people in that run.
NV is the number of random coefficients in that run.

The resulting DR is a 3-dimensional array of standardized draws, where element DR(i,j,k) is the k-th draw for person j for random coefficient i. If you want to use these draws in a new run, then set DRAWTYPE=5 and put the above statements after "if DRAWTYPE==5". Of course, the new run must have the same value of NP, NDRAWS, and NV as the run that created the file, and must set NDRAWS=NTAKES.

If you cannot (or do not want to) hold the entire file in memory at one time, you can read-in only a portion of it. To read K draws of all the random coefficients for all people, use these commands:
m=memmapfile('filename','Format',{'double',[NDRAWS NP NV],'j'});
DR=m.Data(1).j(1:K,:,:);
DR=permute(DR,[3,2,1]);
Then element DR(i,j,k) is the k-th draw for person j for coefficient i.

You might wonder why the permute statement is needed. You actually don't need it if you are just looking at the draws for your own purposes; without the statement, DR is simply ordered differently, such that DR(k,j,i) is the k-th draw for person j for coefficient i. However, if you want to use the draws in with the estimation code, then the permute statement is needed. There is a reason that the draws are stored in a different arrangement than is used in the code, in case you care to know: For Halton draws, an entire sequence must be created together for a given random coefficient, even if NMEM<NDRAWS, since later parts of the sequence are created from earlier parts (see Train, 2003, pp. 225-6.) So, in order to avoid running out of memory when holding all NDRAWSxNP draws, the draws for only *one* random coefficient are created at a time and then written into the file. This process creates a file of numbers with the ordering: all NDRAWSxNP draws for the first coefficient, followed by all NDRAWSxNP draws for the second coefficient, etc. When read in, this file makes a 3-dimensional array with NV as the *last* dimension. However, the code uses draws in an array of dimension NVxNPxNMEM, where NV is the *first* dimension. The permute statement appropriately transposes the first and third dimensions. 




