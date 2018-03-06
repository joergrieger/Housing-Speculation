clear

global exstate endstate1 hh1housetoday xgrid
global FH FL
global gamma1 gamma2 beta transcost margin pix1 pix2
global end1 end2
global hh1bond1 hh1bond2 hh1bond3 hh1bond4 hh1bond5 hh1bond6 hh1bond7 hh1bond8 
global pbond1 pbond2 pbond3 pbond4 pbond5 pbond6 pbond7 pbond8
global hh1house1 hh1house2 hh1house3 hh1house4 hh1house5 hh1house6 hh1house7 hh1house8
global phouse1 phouse2 phouse3 phouse4 phouse5 phouse6 phouse7 phouse8
global Value101 Value102 Value103 Value104 Value105 Value106 Value107 Value108
global Value201 Value202 Value203 Value204 Value205 Value206 Value207 Value208
global r101 r102 r103 r104 r105 r106 r107 r108 
global r201 r202 r203 r204 r205 r206 r207 r208 
global ch101 ch102 ch103 ch104 ch105 ch106 ch107 ch108
global ch201 ch202 ch203 ch204 ch205 ch206 ch207 ch208
global ev1 ev2 psi1 psi2 growth
global alpha
nendog=101;
%---------------------------------------------------------------------------
% Variables for the time iteration
%---------------------------------------------------------------------------
pbondn      =zeros(nendog,nendog,8);
hh1bondn    =zeros(nendog,nendog,8);
hh1housen   =ones(nendog,nendog,8)*0.5;
phousen     =zeros(nendog,nendog,8);
ch1n        =zeros(nendog,nendog,8);
ch2n        =zeros(nendog,nendog,8);
r1n         =zeros(nendog,nendog,8);
r2n         =zeros(nendog,nendog,8);
Value1n     =ones(nendog,nendog,8);
Value2n     =ones(nendog,nendog,8);

pbond       =zeros(nendog,nendog,8);
hh1bond     =zeros(nendog,nendog,8);
hh1house    =ones(nendog,nendog,8)*0.5;
phouse      =ones(nendog,nendog,8)*log(0.05);
ch1         =ones(nendog,nendog,8);
ch2         =ones(nendog,nendog,8);
r1          =ones(nendog,nendog,8);
r2          =ones(nendog,nendog,8);
Value1      =ones(nendog,nendog,8)*0;
Value2      =ones(nendog,nendog,8)*0;
optimSolutions=ones(nendog,nendog,8,24)*(-0.1);

%---------------------------------------------------------------------------
% Parameterize the Economy: Beliefs
%---------------------------------------------------------------------------

a1=0.50;
a2=0.14;
a3=0.14;
a4=0.14;

alpha1=0.57;
alpha2=0.57;
etah=1.6;

A=[a1, alpha1-a1, alpha2-a1, 1+a1-alpha1-alpha2;
   a2, alpha1-a2, alpha2-a2, 1+a2-alpha1-alpha2;
   a3, alpha1-a3, alpha2-a3, 1+a3-alpha1-alpha2;
   a4, alpha1-a4, alpha2-a4, 1+a4-alpha1-alpha2];

transMat=[0.43 0.57;
          0.57 0.43];
Gamma=[transMat(1,1)*A transMat(1,2)*A;
       transMat(2,1)*A transMat(2,2)*A];
   

FH=[   etah*transMat(1,1)*A transMat(1,2)*(1-etah*(transMat(1,1)))/(transMat(1,2))*A;
       etah*transMat(2,1)*A transMat(2,2)*(1-etah*(transMat(2,1)))/(transMat(2,2))*A];

FL=1/(1-alpha1)*(Gamma-alpha1*FH);

%---------------------------------------------------------------------------
% Parameterize the Economy: Preferences
%---------------------------------------------------------------------------
alpha=0.95;
gamma1=2.0;
gamma2=2.0;
psi1=1.5;
psi2=1.5;
beta=0.96;
pix1=1.25;
pix2=1.25;
%---------------------------------------------------------------------------
% Parameterize the Economy: Economic Fundamentals
%---------------------------------------------------------------------------
aggGrowthHigh=1.054;
aggGrowthLow=0.982;
growth=[aggGrowthHigh aggGrowthHigh aggGrowthHigh aggGrowthHigh aggGrowthLow aggGrowthLow aggGrowthLow aggGrowthLow];
margin=0.00;
end1=[0.50 0.50 0.50 0.50 0.50 0.50 0.50 0.50];
end2=1-end1;
transcost=35*1e-4;
%---------------------------------------------------------------------------
% Create grid2*1e-4;

%---------------------------------------------------------------------------
grid=linspace(0.00,1.000,nendog);
hgrid=linspace(0.01,0.99,nendog);
[xgrid,ygrid]=ndgrid(grid,hgrid);

%---------------------------------------------------------------------------
% Set options for solvers
%---------------------------------------------------------------------------
opts=optiset('solver','ipopt','maxiter',10000,'display','iter','tolafun',1e-7,'tolrfun',1e-7);
options = optimoptions(@fmincon,'MaxIter',50000,'MaxFunEvals',5000000,'Display','iter','Algorithm','sqp','SubproblemAlgorithm','cg','FinDiffType','central','TolCon',1e-7);

%---------------------------------------------------------------------------
% Load old solutions
%---------------------------------------------------------------------------
%load('Solutions_tempb.mat');
%load('Solutions_c035_m050_p150_e160.mat');
%--------------------------------------------------------------------------
% Upper and lower bounds
%--------------------------------------------------------------------------
lb=ones(24,1)*(-Inf);
ub=ones(24,1)*(Inf);
lb(11:14,1)=0;
lb(23:24,1)=1-transcost;
ub(23:24,1)=1+transcost;
%---------------------------------------------------------------------------
% Start time iteration
%---------------------------------------------------------------------------
max_error=100;
time_iter=0;
while max_error>1*1e-4
    time_iter=time_iter+1;
    disp(time_iter);
    
    Value101=griddedInterpolant(xgrid,ygrid,Value1(:,:,1));
    Value102=griddedInterpolant(xgrid,ygrid,Value1(:,:,2));
    Value103=griddedInterpolant(xgrid,ygrid,Value1(:,:,3));
    Value104=griddedInterpolant(xgrid,ygrid,Value1(:,:,4));
    Value105=griddedInterpolant(xgrid,ygrid,Value1(:,:,5));
    Value106=griddedInterpolant(xgrid,ygrid,Value1(:,:,6));
    Value107=griddedInterpolant(xgrid,ygrid,Value1(:,:,7));
    Value108=griddedInterpolant(xgrid,ygrid,Value1(:,:,8));
    
    Value201=griddedInterpolant(xgrid,ygrid,Value2(:,:,1));
    Value202=griddedInterpolant(xgrid,ygrid,Value2(:,:,2));
    Value203=griddedInterpolant(xgrid,ygrid,Value2(:,:,3));
    Value204=griddedInterpolant(xgrid,ygrid,Value2(:,:,4));
    Value205=griddedInterpolant(xgrid,ygrid,Value2(:,:,5));
    Value206=griddedInterpolant(xgrid,ygrid,Value2(:,:,6));
    Value207=griddedInterpolant(xgrid,ygrid,Value2(:,:,7));
    Value208=griddedInterpolant(xgrid,ygrid,Value2(:,:,8));
    
    ch101   =griddedInterpolant(xgrid,ygrid,ch1(:,:,1));
    ch102   =griddedInterpolant(xgrid,ygrid,ch1(:,:,2));
    ch103   =griddedInterpolant(xgrid,ygrid,ch1(:,:,3));
    ch104   =griddedInterpolant(xgrid,ygrid,ch1(:,:,4));
    ch105   =griddedInterpolant(xgrid,ygrid,ch1(:,:,5));
    ch106   =griddedInterpolant(xgrid,ygrid,ch1(:,:,6));
    ch107   =griddedInterpolant(xgrid,ygrid,ch1(:,:,7));
    ch108   =griddedInterpolant(xgrid,ygrid,ch1(:,:,8));
    
    ch201   =griddedInterpolant(xgrid,ygrid,ch2(:,:,1));
    ch202   =griddedInterpolant(xgrid,ygrid,ch2(:,:,2));
    ch203   =griddedInterpolant(xgrid,ygrid,ch2(:,:,3));
    ch204   =griddedInterpolant(xgrid,ygrid,ch2(:,:,4));
    ch205   =griddedInterpolant(xgrid,ygrid,ch2(:,:,5));
    ch206   =griddedInterpolant(xgrid,ygrid,ch2(:,:,6));
    ch207   =griddedInterpolant(xgrid,ygrid,ch2(:,:,7));
    ch208   =griddedInterpolant(xgrid,ygrid,ch2(:,:,8));
    
    hh1bond1 = griddedInterpolant(xgrid,ygrid,hh1bond(:,:,1));
    hh1bond2 = griddedInterpolant(xgrid,ygrid,hh1bond(:,:,2));
    hh1bond3 = griddedInterpolant(xgrid,ygrid,hh1bond(:,:,3));
    hh1bond4 = griddedInterpolant(xgrid,ygrid,hh1bond(:,:,4));
    hh1bond5 = griddedInterpolant(xgrid,ygrid,hh1bond(:,:,5));
    hh1bond6 = griddedInterpolant(xgrid,ygrid,hh1bond(:,:,6));
    hh1bond7 = griddedInterpolant(xgrid,ygrid,hh1bond(:,:,7));
    hh1bond8 = griddedInterpolant(xgrid,ygrid,hh1bond(:,:,8));
    
    pbond1 = griddedInterpolant(xgrid,ygrid,pbond(:,:,1));
    pbond2 = griddedInterpolant(xgrid,ygrid,pbond(:,:,2));
    pbond3 = griddedInterpolant(xgrid,ygrid,pbond(:,:,3));
    pbond4 = griddedInterpolant(xgrid,ygrid,pbond(:,:,4));
    pbond5 = griddedInterpolant(xgrid,ygrid,pbond(:,:,5));
    pbond6 = griddedInterpolant(xgrid,ygrid,pbond(:,:,6));
    pbond7 = griddedInterpolant(xgrid,ygrid,pbond(:,:,7));
    pbond8 = griddedInterpolant(xgrid,ygrid,pbond(:,:,8));
    
    hh1house1= griddedInterpolant(xgrid,ygrid,hh1house(:,:,1));
    hh1house2= griddedInterpolant(xgrid,ygrid,hh1house(:,:,2));
    hh1house3= griddedInterpolant(xgrid,ygrid,hh1house(:,:,3));
    hh1house4= griddedInterpolant(xgrid,ygrid,hh1house(:,:,4));
    hh1house5= griddedInterpolant(xgrid,ygrid,hh1house(:,:,5));
    hh1house6= griddedInterpolant(xgrid,ygrid,hh1house(:,:,6));
    hh1house7= griddedInterpolant(xgrid,ygrid,hh1house(:,:,7));
    hh1house8= griddedInterpolant(xgrid,ygrid,hh1house(:,:,8));
    
    phouse1= griddedInterpolant(xgrid,ygrid,phouse(:,:,1));
    phouse2= griddedInterpolant(xgrid,ygrid,phouse(:,:,2));
    phouse3= griddedInterpolant(xgrid,ygrid,phouse(:,:,3));
    phouse4= griddedInterpolant(xgrid,ygrid,phouse(:,:,4));
    phouse5= griddedInterpolant(xgrid,ygrid,phouse(:,:,5));
    phouse6= griddedInterpolant(xgrid,ygrid,phouse(:,:,6));
    phouse7= griddedInterpolant(xgrid,ygrid,phouse(:,:,7));
    phouse8= griddedInterpolant(xgrid,ygrid,phouse(:,:,8));
    
    r101= griddedInterpolant(xgrid,ygrid,r1(:,:,1));
    r102= griddedInterpolant(xgrid,ygrid,r1(:,:,2));
    r103= griddedInterpolant(xgrid,ygrid,r1(:,:,3));
    r104= griddedInterpolant(xgrid,ygrid,r1(:,:,4));
    r105= griddedInterpolant(xgrid,ygrid,r1(:,:,5));
    r106= griddedInterpolant(xgrid,ygrid,r1(:,:,6));
    r107= griddedInterpolant(xgrid,ygrid,r1(:,:,7));
    r108= griddedInterpolant(xgrid,ygrid,r1(:,:,8));
    
    r201= griddedInterpolant(xgrid,ygrid,r2(:,:,1));
    r202= griddedInterpolant(xgrid,ygrid,r2(:,:,2));
    r203= griddedInterpolant(xgrid,ygrid,r2(:,:,3));
    r204= griddedInterpolant(xgrid,ygrid,r2(:,:,4));
    r205= griddedInterpolant(xgrid,ygrid,r2(:,:,5));
    r206= griddedInterpolant(xgrid,ygrid,r2(:,:,6));
    r207= griddedInterpolant(xgrid,ygrid,r2(:,:,7));
    r208= griddedInterpolant(xgrid,ygrid,r2(:,:,8));
    
    %xguess=ones(30,1)*0.501;        
	for exstate=1:8
        disp(exstate);
		xguess=squeeze(optimSolutions(1,1,exstate,:));
        if exstate>1
            xguess=squeeze(optimSolutions(1,1,exstate-1,:));
        end
        housestate=1;
        for housestate=1:nendog
            hh1housetoday=hgrid(housestate);            
            if housestate>1
                xguess=squeeze(optimSolutions(1,housestate-1,exstate,:));
            end
            for endstate1=1:nendog
                if time_iter>1
                    xguess=squeeze(optimSolutions(endstate1,housestate,exstate,:));
                end
                [x,z,exitflag]=opti_fmincon(@FDummy,xguess,[],[],[],[],lb,ub,@EulerEquation,opts);
                xguess=x;
                [x,z,exitflag]=fmincon(@FDummy,xguess,[],[],[],[],[],[],@EulerEquation,options);
                xguess=x;
                %if exitflag<0
                %    pause;
                %end
                %pause;
                pbondn(endstate1,housestate,exstate)=x(21);
                phousen(endstate1,housestate,exstate)=x(22);
                hh1bondn(endstate1,housestate,exstate)=x(9);
                hh1housen(endstate1,housestate,exstate)=hh1housetoday+x(11)-x(12);
                ch1n(endstate1,housestate,exstate)=x(15);
                ch2n(endstate1,housestate,exstate)=x(16);
                r1n(endstate1,housestate,exstate)=x(23);
                r2n(endstate1,housestate,exstate)=x(24);
                Value1n(endstate1,housestate,exstate)=log(ev1);
                Value2n(endstate1,housestate,exstate)=log(ev2);
                optimSolutions(endstate1,housestate,exstate,:)=x;
                %pause;
            end
            %pause;
        end
        save('temporary.mat','exstate','phousen','hh1housen');
        %pause;
    end
    %for ii=1:8
    %    for jj=1:nendog
    %        for kk=2:nendog
    %            pbondn(jj,kk,ii)=pbondn(jj,1,ii);
    %            phousen(jj,kk,ii)=phousen(jj,1,ii);
    %            hh1bondn(jj,kk,ii)=hh1bondn(jj,1,ii);
    %            hh1housen(jj,kk,ii)=hh1housen(jj,1,ii);
    %            ch1n(jj,kk,ii)=ch1n(jj,1,ii);
    %            ch2n(jj,kk,ii)=ch2n(jj,1,ii);
    %            r1n(jj,kk,ii)=r1n(jj,1,ii);
    %            r2n(jj,kk,ii)=r2n(jj,1,ii);
    %            Value1n(jj,kk,ii)=Value1n(jj,1,ii);
    %            Value2n(jj,kk,ii)=Value2n(jj,1,ii);
    %            optimSolutions(jj,kk,ii,:)=squeeze(optimSolutions(jj,1,ii,:));
    %        end
    %    end
    %end
    diff_bond = max(max(max(abs(exp(pbondn)-exp(pbond)))));
    diff_house = max(max(max(abs(exp(phousen)-exp(phouse)))));
    diff_hh1b1= max(max(max(abs((hh1bondn)-(hh1bond)))));
    diff_hh1h1 =max(max(max(abs(hh1housen-hh1house))));
    diff=[diff_bond diff_house diff_hh1b1 diff_hh1h1];
	max_error=max(diff);
	disp('Maximum Differences in Prices and Portfolios');
	disp(diff);
	pbond=pbondn;
    phouse=phousen;
	ch1=ch1n;
	ch2=ch2n;
    Value1=Value1n;
    Value2=Value2n;
	hh1bond=hh1bondn;
    hh1house=hh1housen;
    r1=r1n;
    r2=r2n;

    save('Solutions_tempb.mat','pbond','phouse','hh1bond','hh1house','ch1','ch2','r1','r2','Value1','Value2','optimSolutions','max_error','diff');
    %pause;
    %options = optimset('Display','iter','TolCon',1e-10,'TolX',1e-30,'MaxIter',50000);
end
save('Solutions_c010_m000_p150_e160.mat','pbond','phouse','hh1bond','hh1house','ch1','ch2','r1','r2','Value1','Value2','optimSolutions','max_error');

		