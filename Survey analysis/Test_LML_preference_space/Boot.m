mnhold=zeros((NV+NNR),NReps);
stdhold=zeros(NV,NReps);
covhold=zeros(NV,NV,NReps);

XMAT_Original=XMAT;
PID_Original=XMAT(:,1);
CSID_Original=XMAT(:,2);
clear i r
for xi=1:NReps
    btsample=randi(NP,[NP,1]);
    XMAT=[];
        for r=1:NP
              thisp=btsample(r,1);
              thisx=XMAT_Original(PID_Original==thisp,:);
              thisx(:,1)=r.*ones(size(thisx,1),1);
              XMAT=[XMAT ; thisx];   
        end
        XMAT(:,2)=CSID_Original;
    DrawsBeta
    CreateZ
    param=StartB;
    options=optimset('LargeScale','off','Display','iter','GradObj','on',...
    'MaxFunEvals',10000,'MaxIter',2000,'TolX',10^(-6),'TolFun',10^(-6),'DerivativeCheck','off');
    [paramhat,fval,exitflag]=fminunc(@flexll,param,options);
	fvalhold(xi) = fval;
    paramhold(:,xi)=paramhat;
    [mnhold(:,xi),stdhold(:,xi),covhold(:,:,xi),freqhold(:,:,xi),midhold(:,:,xi)]=stats(paramhat,NBins);    
end
MeanSE=std(mnhold,0,2);
StdSE=std(stdhold,0,2);
Cov_SE=std(covhold,0,3); 

disp('Utility Parameters');
disp(' ');
disp('                    Mean           ');
disp('              ------------------  ');
disp('                 Est     SE      ');
for r=1:length(NAMES)
    fprintf('%-10s %10.4f %10.4f\n', NAMES{r,1}, [MeanEst(r,1) MeanSE(r,1)]);
end

disp(' ');
disp('Estimated var-covar matrix');
disp(CovMatEst);
disp('std. error of var-covar');
disp(Cov_SE);


