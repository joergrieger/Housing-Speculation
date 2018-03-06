function [fvalin,fval] = EulerEquation(x)
% EulerEquation
% Calculates the first order conditions

global exstate endstate1 hh1housetoday xgrid
global FH FL
global gamma1 gamma2 beta transcost margin pix1 pix2
global end1
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

fval=zeros(24,1);
fvalin=zeros(24,1);
%---------------------------------------------------------
% first set of variables: guesses of future wealth share
%---------------------------------------------------------
hh1wshareguess=zeros(8,1);
hh1wshareguess(1:8)=x(1:8);

%-------------------------------------------------------
% guesses for portfolio choices
%-------------------------------------------------------
hh1bondguess        =x(9);
hh2bondguess        =x(10);
hh1housebuyguess    =x(11);
hh1housesellguess   =x(12);
hh1houseguess       =hh1housetoday+hh1housebuyguess-hh1housesellguess;
hh2housebuyguess    =x(13);
hh2housesellguess   =x(14);
hh2houseguess       =(1-hh1housetoday)+hh2housebuyguess-hh2housesellguess;

%------------------------------------------------------
% guesses for consumption
%------------------------------------------------------
ch1=exp(x(15));
ch2=exp(x(16));
%------------------------------------------------------
% guesses for lagrange multiplier
%------------------------------------------------------
lag1collpos         =max(0,x(17))^2;
lag1collneg         =max(0,-x(17))^2;
lag2collpos         =max(0,x(18))^2;
lag2collneg         =max(0,-x(18))^2;

lag1houseshortpos  =max(0,x(19))^2;
lag1houseshortneg  =max(0,-x(19))^2;
lag2houseshortpos  =max(0,x(20))^2;
lag2houseshortneg  =max(0,-x(20))^2;

%-------------------------------------------------------
% guesses for prices
%------------------------------------------------------
qbg=exp(x(21));
qhg=exp(x(22));
%------------------------------------------------------
% guesses for r
%------------------------------------------------------
r1=x(23);
r2=x(24);
%---------------------------------------------------------
% Interpolate future prices and consumption
%---------------------------------------------------------
hh1houset=zeros(8,1);
phouset=zeros(8,1);
Value1t=zeros(8,1);
Value2t=zeros(8,1);
ch1t=zeros(8,1);
ch2t=zeros(8,1);

ii=1;
hh1houset(ii)       =hh1house1(hh1wshareguess(ii),hh1houseguess);
phouset(ii)         =phouse1(hh1wshareguess(ii),hh1houseguess);
Value1t(ii)         =Value101(hh1wshareguess(ii),hh1houseguess);
Value2t(ii)         =Value201(hh1wshareguess(ii),hh1houseguess);
ch1t(ii)            =exp(ch101(hh1wshareguess(ii),hh1houseguess));
ch2t(ii)            =exp(ch201(hh1wshareguess(ii),hh1houseguess));
r1t(ii)             =r101(hh1wshareguess(ii),hh1houseguess);
r2t(ii)             =r201(hh1wshareguess(ii),hh1houseguess);

ii=2;
hh1houset(ii)       =hh1house2(hh1wshareguess(ii),hh1houseguess);
phouset(ii)         =phouse2(hh1wshareguess(ii),hh1houseguess);
Value1t(ii)         =Value102(hh1wshareguess(ii),hh1houseguess);
Value2t(ii)         =Value202(hh1wshareguess(ii),hh1houseguess);
ch1t(ii)            =exp(ch102(hh1wshareguess(ii),hh1houseguess));
ch2t(ii)            =exp(ch202(hh1wshareguess(ii),hh1houseguess));
r1t(ii)             =r102(hh1wshareguess(ii),hh1houseguess);
r2t(ii)             =r202(hh1wshareguess(ii),hh1houseguess);

ii=3;
hh1houset(ii)       =hh1house3(hh1wshareguess(ii),hh1houseguess);
phouset(ii)         =phouse3(hh1wshareguess(ii),hh1houseguess);
Value1t(ii)         =Value103(hh1wshareguess(ii),hh1houseguess);
Value2t(ii)         =Value203(hh1wshareguess(ii),hh1houseguess);
ch1t(ii)            =exp(ch103(hh1wshareguess(ii),hh1houseguess));
ch2t(ii)            =exp(ch203(hh1wshareguess(ii),hh1houseguess));
r1t(ii)             =r103(hh1wshareguess(ii),hh1houseguess);
r2t(ii)             =r203(hh1wshareguess(ii),hh1houseguess);

ii=4;
hh1houset(ii)       =hh1house4(hh1wshareguess(ii),hh1houseguess);
phouset(ii)         =phouse4(hh1wshareguess(ii),hh1houseguess);
Value1t(ii)         =Value104(hh1wshareguess(ii),hh1houseguess);
Value2t(ii)         =Value204(hh1wshareguess(ii),hh1houseguess);
ch1t(ii)            =exp(ch104(hh1wshareguess(ii),hh1houseguess));
ch2t(ii)            =exp(ch204(hh1wshareguess(ii),hh1houseguess));
r1t(ii)             =r104(hh1wshareguess(ii),hh1houseguess);
r2t(ii)             =r204(hh1wshareguess(ii),hh1houseguess);

ii=5;
hh1houset(ii)       =hh1house5(hh1wshareguess(ii),hh1houseguess);
phouset(ii)         =phouse5(hh1wshareguess(ii),hh1houseguess);
Value1t(ii)         =Value105(hh1wshareguess(ii),hh1houseguess);
Value2t(ii)         =Value205(hh1wshareguess(ii),hh1houseguess);
ch1t(ii)            =exp(ch105(hh1wshareguess(ii),hh1houseguess));
ch2t(ii)            =exp(ch205(hh1wshareguess(ii),hh1houseguess));
r1t(ii)             =r105(hh1wshareguess(ii),hh1houseguess);
r2t(ii)             =r205(hh1wshareguess(ii),hh1houseguess);

ii=6;
hh1houset(ii)       =hh1house6(hh1wshareguess(ii),hh1houseguess);
phouset(ii)         =phouse6(hh1wshareguess(ii),hh1houseguess);
Value1t(ii)         =Value106(hh1wshareguess(ii),hh1houseguess);
Value2t(ii)         =Value206(hh1wshareguess(ii),hh1houseguess);
ch1t(ii)            =exp(ch106(hh1wshareguess(ii),hh1houseguess));
ch2t(ii)            =exp(ch206(hh1wshareguess(ii),hh1houseguess));
r1t(ii)             =r106(hh1wshareguess(ii),hh1houseguess);
r2t(ii)             =r206(hh1wshareguess(ii),hh1houseguess);

ii=7;
hh1houset(ii)       =hh1house7(hh1wshareguess(ii),hh1houseguess);
phouset(ii)         =phouse7(hh1wshareguess(ii),hh1houseguess);
Value1t(ii)         =Value107(hh1wshareguess(ii),hh1houseguess);
Value2t(ii)         =Value207(hh1wshareguess(ii),hh1houseguess);
ch1t(ii)            =exp(ch107(hh1wshareguess(ii),hh1houseguess));
ch2t(ii)            =exp(ch207(hh1wshareguess(ii),hh1houseguess));
r1t(ii)             =r107(hh1wshareguess(ii),hh1houseguess);
r2t(ii)             =r207(hh1wshareguess(ii),hh1houseguess);

ii=8;
hh1houset(ii)       =hh1house8(hh1wshareguess(ii),hh1houseguess);
phouset(ii)         =phouse8(hh1wshareguess(ii),hh1houseguess);
Value1t(ii)         =Value108(hh1wshareguess(ii),hh1houseguess);
Value2t(ii)         =Value208(hh1wshareguess(ii),hh1houseguess);
ch1t(ii)            =exp(ch108(hh1wshareguess(ii),hh1houseguess));
ch2t(ii)            =exp(ch208(hh1wshareguess(ii),hh1houseguess));
r1t(ii)             =r108(hh1wshareguess(ii),hh1houseguess);
r2t(ii)             =r208(hh1wshareguess(ii),hh1houseguess);

Value1t=exp(Value1t);
Value2t=exp(Value2t);
phouset=exp(phouset);
pminhouse=min(phouset);
hh2houset=1-hh1houset;

ch1tilde=zeros(8,1);
ch2tilde=zeros(8,1);
for ii=1:8
    ch1tilde(ii) = (alpha*ch1t(ii)^((pix1-1)/pix1)+(1-alpha)*hh1houset(ii)^((pix1-1)/pix1))^(pix1/(pix1-1));
    ch2tilde(ii) = (alpha*ch2t(ii)^((pix2-1)/pix2)+(1-alpha)*hh2houset(ii)^((pix2-1)/pix2))^(pix2/(pix2-1));
end
%--------------------------------------------------------
% First order conditions
%--------------------------------------------------------

qh1=0;
qh2=0;
qb1=0;
qb2=0;
switch exstate
    case 1
        transMatrix1=FH;
        transMatrix2=FH;
    case 2
        transMatrix1=FH;
        transMatrix2=FL;
    case 3
        transMatrix1=FL;
        transMatrix2=FH;
    case 4
        transMatrix1=FL;
        transMatrix2=FL;
    case 5
        transMatrix1=FH;
        transMatrix2=FH;
    case 6
        transMatrix1=FH;
        transMatrix2=FL;
    case 7
        transMatrix1=FL;
        transMatrix2=FH;
    case 8
        transMatrix1=FL;
        transMatrix2=FL;
end
EValue1=0;
EValue2=0;
for ii=1:8
    qb1 = qb1+transMatrix1(exstate,ii)*Value1t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma1)*(1-beta)*ch1tilde(ii)^(-1/psi1)*ch1tilde(ii)^(1/pix1)*alpha*ch1t(ii)^(1/(1-pix1))/growth(ii);
    qb2 = qb2+transMatrix2(exstate,ii)*Value2t(ii)^(1/psi2-gamma1)*growth(ii)^(1-gamma2)*(1-beta)*ch2tilde(ii)^(-1/psi2)*ch2tilde(ii)^(1/pix2)*alpha*ch2t(ii)^(1/(1-pix2))/growth(ii);

    qh1=qh1+transMatrix1(exstate,ii)*Value1t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma1)*(1-beta)*ch1tilde(ii)^(-1/psi2)*ch1tilde(ii)^(1/pix1)*alpha*ch1t(ii)^(1/(1-pix1))*(r1t(ii)*phouset(ii));
    qh2=qh2+transMatrix2(exstate,ii)*Value2t(ii)^(1/psi1-gamma1)*growth(ii)^(1-gamma2)*(1-beta)*ch2tilde(ii)^(-1/psi2)*ch2tilde(ii)^(1/pix2)*alpha*ch2t(ii)^(1/(1-pix2))*(r2t(ii)*phouset(ii));

    EValue1=EValue1+transMatrix1(exstate,ii)*(Value1t(ii)*growth(ii))^(1-gamma1);
    EValue2=EValue2+transMatrix2(exstate,ii)*(Value2t(ii)*growth(ii))^(1-gamma2);
end
ch1til=(alpha*ch1^((pix1-1)/pix1)+(1-alpha)*hh1houseguess^((pix1-1)/pix1))^(pix1/(pix1-1));
ch2til=(alpha*ch2^((pix2-1)/pix2)+(1-alpha)*hh2houseguess^((pix2-1)/pix2))^(pix2/(pix2-1));

phi1=(1-gamma1)/(1-1/psi1);
phi2=(1-gamma2)/(1-1/psi2);

Value1today=((1-beta)*ch1til^(1-1/psi1)+beta*EValue1^(1/phi1))^(phi1/(1-gamma1));
Value2today=((1-beta)*ch2til^(1-1/psi2)+beta*EValue2^(1/phi2))^(phi2/(1-gamma2));

MValue1today=(1-beta)*Value1today^(1/psi1)*ch1til^(-1/psi1)*alpha*ch1til^(1/pix1)*ch1^(1/(1-pix1));
MValue2today=(1-beta)*Value2today^(1/psi2)*ch2til^(-1/psi2)*alpha*ch2til^(1/pix2)*ch2^(1/(1-pix2));

MValueh1today=(1-beta)*Value1today^(1/psi1)*ch1til^(-1/psi1)*(1-alpha)*ch1til^(1/pix1)*hh1houseguess^(1/(1-pix1));
MValueh2today=(1-beta)*Value2today^(1/psi2)*ch2til^(-1/psi2)*(1-alpha)*ch2til^(1/pix1)*hh2houseguess^(1/(1-pix2));

ev1=Value1today;
ev2=Value2today;
%--------------------------------------------------------
% Calculate financial wealth share
%-------------------------------------------------------
pw=hh1houseguess*r1t+(1-hh1houseguess)*r2t;
    hh1wsharet=(hh1houseguess*phouset.*r1t'+hh1bondguess.*1./growth')./(phouset.*pw');
    fval(1:8)=hh1wsharet-hh1wshareguess;
%---------------------------------------------------------------------
% First order conditions
%---------------------------------------------------------------------
fval(9)=qhg*r1*MValue1today-(lag1houseshortpos+lag1collpos*margin*pminhouse+MValueh1today+beta*Value1today^(1/psi1)*EValue1^((1-phi1)/phi1)*qh1);
fval(10)=qhg*r2*MValue2today-(lag2houseshortpos+lag2collpos*margin*pminhouse+MValueh2today+beta*Value2today^(1/psi2)*EValue2^((1-phi2)/phi2)*qh2);

fval(11)=qbg*MValue1today-(lag1collpos+beta*Value1today^(1/psi1)*EValue1^((1-phi1)/phi1)*qb1);
fval(12)=qbg*MValue2today-(lag2collpos+beta*Value2today^(1/psi2)*EValue2^((1-phi2)/phi2)*qb2);
%----------------------------------------------------------------------
% Budget constraint
%----------------------------------------------------------------------
priceweight=hh1housetoday*r1+(1-hh1housetoday)*r2;
taxincome=(hh1housetoday-hh1houseguess)*r1*qhg+((1-hh1housetoday)-hh1houseguess)*r2*qhg;
income1=end1(exstate)+xgrid(endstate1)*(priceweight*qhg)-hh1bondguess*qbg-hh1houseguess*r1*qhg+0.5*taxincome;
fval(13)=ch1-income1;
fval(14)=1-ch1-ch2;
fval(15)=1-hh1houseguess-hh2houseguess;
fval(16)=hh1bondguess+hh2bondguess;
%-------------------------------------------------------------------
% Complementary slackness conditions
%-------------------------------------------------------------------
fval(17)=lag1collneg-(hh1houseguess*pminhouse*margin+hh1bondguess);
fval(18)=lag2collneg-(hh2houseguess*pminhouse*margin+hh2bondguess);
fval(19)=(1+transcost-r1)*(hh1housesellguess);
fval(20)=(1+transcost-r2)*(hh2housesellguess);
fval(21)=(1-transcost-r1)*(hh1housebuyguess);
fval(22)=(1-transcost-r2)*(hh2housebuyguess);
fval(23)=lag1houseshortneg-hh1houseguess;
fval(24)=lag2houseshortneg-hh2houseguess;
%------------------------------------------------------------------
% Inequality constraints
%------------------------------------------------------------------
fvalin(1)=1-transcost-r1;
fvalin(2)=1-transcost-r2;
fvalin(3)=r1-(1+transcost);
fvalin(4)=r2-(1+transcost);
fvalin(5)=-hh1housebuyguess;
fvalin(6)=-hh2housebuyguess;
fvalin(7)=-hh1housesellguess;
fvalin(8)=-hh2housesellguess;
end
