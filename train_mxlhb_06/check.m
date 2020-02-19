% This code checks the input data and specifications that the user provides in mxlhb.m
% Written by Kenneth Train, July 26, 2006, latest edits on Sept 24, 2006.

function ok=check

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

% Check for positive intergers

ok=0;

if ceil(NP) ~= NP | NP < 1;
   disp(['NP must be a positive integer, but it is set to ' num2str(NP)]);
   disp('Program terminated.');
   return
end

if ceil(NCS) ~= NCS | NCS < 1;
   disp(['NCS must be a positive integer, but it is set to ' num2str(NCS)]);
   disp('Program terminated.');
   return
end

if ceil(NROWS) ~= NROWS | NROWS < 1;
   disp(['NROWS must be a positive integer, but it is set to ' num2str(NROWS)]);
   disp('Program terminated.');
   return
end

if ceil(NCREP) ~= NCREP | NCREP < 1;
   disp(['NCREP must be a positive integer, but it is set to ' num2str(NCREP)]);
   disp('Program terminated.');
   return
end

if ceil(NEREP) ~= NEREP | NEREP < 1;
   disp(['NEREP must be a positive integer, but it is set to ' num2str(NEREP)]);
   disp('Program terminated.');
   return
end

if ceil(NSKIP) ~= NSKIP | NSKIP < 1;
   disp(['NSKIP must be a positive integer, but it is set to ' num2str(NSKIP)]);
   disp('Program terminated.');
   return
end

if ceil(NRFOR) ~= NRFOR | NRFOR < 1;
   disp(['NRFOR must be a positive integer, but it is set to ' num2str(NRFOR)]);
   disp('Program terminated.');
   return
end

if ceil(NRLL) ~= NRLL | NRLL < 1;
   disp(['NRLL must be a positive integer, but it is set to ' num2str(NRLL)]);
   disp('Program terminated.');
   return
end

if ceil(INFOSKIP) ~= INFOSKIP | INFOSKIP < 1;
   disp(['INFOSKIP must be a positive integer, but it is set to ' num2str(INFOSKIP)]);
   disp('Program terminated.');
   return
end

if ceil(SEED1) ~= SEED1 | SEED1 < 1;
   disp(['SEED1 must be a positive integer, but it is set to ' num2str(SEED1)]);
   disp('Program terminated.');
   return
end

if ceil(SEED2) ~= SEED2 | SEED2 < 1;
   disp(['SEED2 must be a positive integer, but it is set to ' num2str(SEED2)]);
   disp('Program terminated.');
   return
end

if NV>0 & sum(sum(ceil(IDV) ~= IDV | IDV < 1),2) ~ 0;
   disp('IDV must contain positive integers only, but it contains other values.');
   disp('Program terminated.');
   return
end

if NF>0 & sum(sum(ceil(IDF) ~= IDF | IDF < 1),2) ~ 0;
   disp('IDF must contain positive integers only, but it contains other values.');
   disp('Program terminated.');
   return
end

if FULLCV ~= 1 & FULLCV ~= 0;
   disp(['FULLCV must be 0 or 1, but it is set to ' num2str(FULLCV)]);
   disp('Program terminated.');
   return
end

if KEEPMN ~= 1 & KEEPMN ~= 0;
   disp(['KEEPMN must be 0 or 1, but it is set to ' num2str(KEEPMN)]);
   disp('Program terminated.');
   return
end

if WANTINDC ~= 1 & WANTINDC ~= 0;
   disp(['WANTINDC must be 0 or 1, but it is set to ' num2str(WANTINDC)]);
   disp('Program terminated.');
   return
end

% Checking Rhos.

if RHO <= 0 ;
   disp(['RHO must be strictly positive, but it is set to ' num2str(RHO)]);
   disp('Program terminated.');
   return
end

if RHOF <= 0 ;
   disp(['RHOF must be strictly positive, but it is set to ' num2str(RHOF)]);
   disp('Program terminated.');
   return
end


% Checking XMAT %
if ( size(XMAT,1) ~= NROWS)
      disp(['XMAT has ' num2str(size(XMAT,1)) ' rows']);
      disp(['but it should have NROWS= '  num2str(NROWS)   ' rows.']);
      disp('Program terminated.');
      return
end

if sum(XMAT(:,1) > NP) ~= 0
     disp(['The first column of XMAT has a value greater than NP= ' num2str(NP)]);
     disp('Program terminated.');
     return
end

if sum(XMAT(:,1) < 1) ~= 0
     disp('The first column of XMAT has a value less than 1.');
     disp('Program terminated.');
     return
end

k=(XMAT(2:NROWS,1) ~= XMAT(1:NROWS-1,1)) & (XMAT(2:NROWS,1) ~= (XMAT(1:NROWS-1,1)+1));
if sum(k) ~= 0
    disp('The first column of XMAT does not ascend from 1 to NP.');
    disp('Program terminated.')
    return
end

if sum(XMAT(:,2) > NCS) ~= 0
     disp(['The second column of XMAT has a value greater than NCS= ' num2str(NCS)]);
     disp('Program terminated.');
     return
end

if sum(XMAT(:,2) < 1) ~= 0
     disp('The second column of XMAT has a value less than 1.');
     disp('Program terminated.');
     return
end

k=(XMAT(2:NROWS,2) ~= XMAT(1:NROWS-1,2)) & (XMAT(2:NROWS,2) ~= (XMAT(1:NROWS-1,2)+1));
if sum(k) ~= 0
    disp('The second column of XMAT does not ascend from 1 to NCS.');
    disp('Program terminated.')
    return
end

if sum(XMAT(:,3) ~= 0 & XMAT(:,3) ~= 1) ~= 0
     disp('The third column of XMAT has a value other than 1 or 0.');
     disp('Program terminated.');
     return
end

for s=1:NCS
    k=(XMAT(:,2) == s);
    if sum(XMAT(k,3)) > 1
       disp('The third column of XMAT indicates more than one chosen alternative');
       disp(['for choice situation ' num2str(s)]);
       disp('Program terminated.');
       return
    end
    if sum(XMAT(k,3)) < 1
       disp('The third column of XMAT indicates that no alternative was chosen');
       disp(['for choice situation ' num2str(s)]);
       disp('Program terminated.');
       return
    end 
end

if sum(sum(isnan(XMAT)),2) ~= 0
   disp('XMAT contains missing data.');
   disp('Program terminated.');
   return
end;

if sum(sum(isinf(XMAT)),2) ~= 0
   disp('XMAT contains an infinite value.');
   disp('Program terminated.');
   return
end;

if NV>0 & size(IDV,2) ~= 2;
   disp(['IDV must have 2 columns and yet it is set to have ' num2str(size(IDV,2))]);
   disp('Program terminated.');
   return
end;

if NV>0 & sum(IDV(:,1) > size(XMAT,2)) ~= 0;
   disp('IDV identifies a variable that is outside XMAT.');
   disp('The first column of IDV is');
   IDV(:,1)
   disp('when each element of this column must be no greater than')
   disp([num2str(size(XMAT,2)) ' which is the number of columns in XMAT.']);
   disp('Program terminated.');
   return
end;

if NV>0 & sum(IDV(:,1) <= 3) ~= 0;
   disp('Each element in the first column of IDV must exceed 3');
   disp('since the first three variables in XMAT cannot be explanatory variables.');
   disp('But the first column of IDV is');
   IDV(:,1)
   disp('which has an element below 3.')
   disp('Program terminated.');
   return
end;

if NV>0 & sum(IDV(:,2) < 0 | IDV(:,2) > 5) ~= 0;
   disp('The second column of IDV must be integers 1-5 identifying the distributions.');
   disp('But the second column of IDV is specified as');
   IDV(:,2)
   disp('which contains a number other than 1-5.')
   disp('Program terminated.');
   return
end;


if NV>0 & size(NAMES,2) ~= 1;
   disp(['NAMES must have 1 columns and yet it is set to have ' num2str(size(NAMES,2))]);
   disp('Be sure to separate names by semicolons.');
   disp('Program terminated.');
   return
end;

if NV>0 & size(IDV,1) ~= size(NAMES,1);
   disp(['IDV and NAMES must have the same length but IDV has length ' num2str(size(IDV,1))]);
   disp(['while NAMES has length ' num2str(size(NAMES,1))]);
   disp('Program terminated.');
   return
end; 

if NV>0 & size(A,2) ~= 1;
   disp(['A must have 1 column and yet it is set to have ' num2str(size(A,2))]);
   disp('Be sure to separate values by semicolons.');
   disp('Program terminated.');
   return
end;
  
if NV>0 & size(A,1) ~= size(IDV,1);
   disp(['A must have the same length as IDV but instead has length ' num2str(size(A,1))]);
   disp('Program terminated.');
   return
end; 

if NV>0 & size(B,1) ~= size(IDV,1);
   disp(['B must have the same length as IDV but instead has length ' num2str(size(B,1))]);
   disp('Program terminated.');
   return
end;

if NV>0 & size(B,2) ~= NP;
   disp(['B must have NP columns but instead has length ' num2str(size(B,2)) ' columns']);
   disp(['while NP is ' num2str(NP)]);
   disp('Program terminated.');
   return
end;

if NV>0 & size(D,1) ~= size(IDV,1) | size(D,2) ~= size(IDV,1);
   disp('D must have as many columns and rows as there are explanatory variables in IDV.');
   disp(['But D has ' num2str(size(D,1)) ' rows and ' num2str(size(D,2)) ' columns,']);
   disp(['while IDV specifies ' num2str(size(IDV,1)) ' variables.']);
   disp('Program terminated.');
   return
end;

if NF>0 & size(IDF,2) ~= 1;
   disp(['IDF must have 1 column and yet it is set to have ' num2str(size(IDF,2))]);
   disp('Be sure to separate elements by semicolons.');
   disp('Program terminated.');
   return
end;

if NF>0 & sum(IDF > size(XMAT,2)) ~= 0;
   disp('IDF identifies a variable that is outside XMAT.');
   disp('IDF is');
   IDF
   disp('when each element must be no greater than');
   disp([num2str(size(XMAT,2)) ' which is the number of columns in XMAT.']);
   disp('Program terminated.');
   return
end;

if NF>0 & sum(IDF <= 3) ~= 0;
   disp('Each element of IDF must exceed 3 since the first three variables');
   disp('of XMAT cannot be explanatory variables.');
   disp('But IDV is');
   IDV(:,1)
   disp('which contains an element below 3.')
   disp('Program terminated.');
   return
end;

if NF>0 & size(NAMESF,2) ~= 1;
   disp(['NAMESF must have 1 column and yet it is set to have ' num2str(size(NAMESF,2))]);
   disp('Be sure to separate names by semicolons.');
   disp('Program terminated.');
   return
end;

if NF>0 & size(IDF,1) ~= size(NAMESF,1);
   disp(['IDF and NAMESF must have the same length but IDF has length ' num2str(size(IDF,1))]);
   disp(['while NAMESF has length ' num2str(size(NAMESF,1))]);
   disp('Program terminated.');
   return
end; 

if NF>0 & size(F,2) ~= 1;
   disp(['F must have 1 column and yet it is set to have ' num2str(size(F,2))]);
   disp('Be sure to separate values by semicolons.');
   disp('Program terminated.');
   return
end;
  
if NF>0 & size(F,1) ~= size(IDF,1);
   disp(['F must have the same length as IDF but instead has length ' num2str(size(F,1))]);
   disp('Program terminated.');
   return
end; 

ok=1;