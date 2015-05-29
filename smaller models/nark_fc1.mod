%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  New Keynesian model with Unit root shock process  %%
%%  2014.4.16                                        %%
%%  Written by Byoungho Bae                 
%%  variable names are modified 2015.2.6 by Yi          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var 
    c, r, w, n, mc, k, y, i, rk, Zmc, dZtr, Zw,
    Ystar, Rstar, PIstar, q_f, rf, f, Zcp
;

varexo
       er, 
       ea, 
       etr, 
       ei, 
       ec, 
       epi, 
       enw, 
       epif, 
       ex, 
       eg, 
       epr, 
       eystar, 
       erstar, 
       epistar, 
       ecp, 
       ew,
       emc,
       epx; 

////////////// define parameters //////////////
parameters  BETA, DELTA, CHIn, ALPHAk,  CHIc, PSIh, PSIf, THETAh, THETAf, TAYLORpi, 
            TAYLORy, TAYLORr, UIPy, UIPf, PHIx, XIh, KAPPAi, GAMMAtr, SP, KAPPArp, 
            ALPHAc, XIc, ALPHAi, XIi, ETA, IOTA, IT, THETAw, 
            PSIw, v, GAMMA1, GAMMA2, UIPr, ALPHAx, XIx, RHOtr, RHOa, RHOi, 
            RHOc,  RHOpi, RHOnw, RHOpif, RHOx, RHOpii, Gss, RHOg,
            RHOw, RHOcp, RHOmc, SIGMAr, SIGMAa, SIGMAi, SIGMAc,  SIGMApi,
            SIGMAnw, SIGMApif, SIGMAx, SIGMAg, SIGMAw,  SIGMAcp,  
            a11, a12, a13, a21, a22, a23, a31, a32, a33, PHIg;

BETA     = 0.999;                // 할인률
DELTA     = 0.025;                // 자본 감가상각률
GAMMAtr   = 0.0095;               // 영구적 기술충격내 성장률
CHIc    = 0.6909;               // 소비습관(consumption habit) 계수
KAPPAi    = 3.339;                // 투자조정 계수
CHIn     = 2;                    // parameter on labor function
ALPHAk    = 0.33;                 // 생산함수내 자본투여 비율
THETAh     = 0.807;                // 국내재 가격 조정 확률 (1- theta)
THETAf    = 0.546;                // 수입재국내재 가격 조정 확률 (1- theta)
THETAw    = 0.541;                // 임금 조정확률
TAYLORpi    = 0.4019;               // 테일러 준칙내 인플레이션 반응 계수
TAYLORy    = 0.312;                // 테일러 준칙내 생산갭 반응 계수
TAYLORr    = 0.9407;               // 테일러 준칙내 금리 지속성 계수
PSIh     = 6.454;                // 국내재 대체탄력성
PSIf    = 5.175;                // 수입재 대체탄력성
PSIw    = 6.091;                // 가계 노동공급의 대체탄력성   
IOTA    = 0.3162;               // Price indexation 계수 
PHIx   = 0.4880;               // 수출함수 지속성 계수
XIh  = 0.3945;               // 수출함수내 환율 탄력성 계수

UIPy     = 0.165;                 // UIP 조건내 계수
UIPf     = 0.03;                 // UIP 조건내 해외부채에 대한 탄력성 계수
UIPr    = 0.5;
ETA     = 1.2;                  // 노동함수내 Frisch 계수
SP  = 0.975;                // 기업가 생존확률
KAPPArp   = 0.05;                 // 리스크 프리미엄 반응계수
ALPHAc   = 0.56723;               // 소비복합재 국내재 비율
XIc   = 3.955;                // 복합소비재 대체 탄력성
ALPHAi   = 0.35088;               // 투자복합재 국내재 비율 
XIi   = 1.673;                // 투자복합재 대체 탄력성
Gss    = 0.20; 
IT = 1.005;              // 물가목표수준 (1.005^4 = 1.020)   
v       = 1;                    // 운전자금(working capital) 대출 비율
GAMMA1     = 0.02;                 // capital utilization 계수
GAMMA2    = 1.2;                  // capital utilization 계수
ALPHAx    = 0.45550;
XIx    = 1.488;
PHIg   = 1;

RHOtr  = 0.8180;                  // 충격구조 지속성 계수
RHOa    = 0.9459;
RHOi  = 0.9273;
RHOc    = 0.6902;
RHOpi  = 0.3224;
RHOnw    = 0.7981;
RHOpif = 0.9325;
RHOx    = 0.8402;
RHOpii = 0.6705;
RHOg    = 0.8052;
RHOw    = 0.4684;
RHOcp   = 0.90716;
RHOmc   = 0.0;

SIGMAr   = 0.0014;             // 충격 표준편차
SIGMAa   = 0.0091;    
SIGMAi = 0.1186;
SIGMAc   = 0.0261;
SIGMApi = 0.0122;
SIGMAnw   = 0.0034;
SIGMApif= 0.0502;
SIGMAx   = 0.0291;
SIGMAg   = 0.0077;
SIGMAw   = 0.188501;
SIGMAcp  = 0.0088;
                            
a11 =  0.8098;                // Foreign VAR 계수  
a12 =  0.0226;
a13 =  -0.0408;
a21 =  0.0374;
a22 =  0.9396;
a23 = -0.0176;
a31 =  0.1177;
a32 =  0.0294;
a33 =  0.5523;



/////////////////////////////////////////////
//               Model (F.O.C)             //
/////////////////////////////////////////////

model;

/////////////////////////////////// shock processes //////////////////////
// permanent tech. shock: ln(A_t) = ln(A_t-1) + GAMMAtr*exp(Zg) + ea //
dZtr    = RHOtr *dZtr(-1)   + (1-RHOtr)*GAMMAtr + etr;                                                
Zmc     = (1-RHOmc)*PSIh    + RHOmc   *Zmc(-1)   + emc;
Zw      = RHOw  *Zw(-1)     + ew;
Zcp     = RHOcp *Zcp(-1)    + ecp;


// fcnb   ////// Flexible-Closed-No FA block: c, r, w, k, n, mc, y, rk, i

exp(dZtr(+1))/(exp(c)-CHIc*exp(c(-1))/exp(dZtr))
= BETA/(exp(c(+1))-CHIc*exp(c)/exp(dZtr(+1)))*(exp(r));

exp(w) = PSIw/(PSIw-1)*(exp(c)-CHIc*exp(c(-1))/exp(dZtr))*exp(n)^ETA;

1/exp(mc) = Zmc/(Zmc-1);

exp(rk) = exp(r) - (1-DELTA);

exp(k) - (1-DELTA)*exp(k(-1))/exp(dZtr) = exp(i);

exp(w) = exp(rk)*exp(Zw)*(1-ALPHAk)/ALPHAk*exp(k(-1))/exp(n)/exp(dZtr);

exp(mc) = (exp(w)/(1-ALPHAk))^(1-ALPHAk)*(exp(rk)/ALPHAk)^ALPHAk;

exp(y) = exp(k(-1))^ALPHAk*exp(n)^(1-ALPHAk)*exp(dZtr)^(-ALPHAk);

exp(y) = exp(c) + exp(i);

(exp(rk) + (1-DELTA))*exp(PIstar)*exp(q_f) = exp(rf)*exp(q_f(+1));
exp(rf) = exp(- UIPf*(exp(q_f)*exp(f)-(exp(steady_state(y))*UIPy))
              - UIPr*(exp(Rstar)-exp(steady_state(Rstar))
                      - (exp(r)-exp(steady_state(r)))))
          *exp(Zcp)*(exp(Rstar));
exp(q_f)*exp(rf(-1))*exp(f(-1)) = exp(q_f)*exp(f)*exp(dZtr)*exp(PIstar);


////fcne

////////////////////////////// Foreign VAR/////////////////////////////////
Ystar   = a11*Ystar(-1)  + a12*Rstar(-1) + a13*PIstar(-1) + eystar;
Rstar   = a21*Ystar(-1)  + a22*Rstar(-1) + a23*PIstar(-1) + erstar;
PIstar  = a31*Ystar(-1)  + a32*Rstar(-1) + a33*PIstar(-1) + epistar;

end;

initval;
r = 0.7056;
c=      		 1.02618;
w =     		 0.44562;
n  =    		 0.327734;
mc  =   		 -0.168351;
k    =  		 3.41132;
y     = 		 1.34218;
i      =		 0.036335;
rk     	=	 -3.33665;
Zmc    	=	 6.454;
dZtr   	=	 0.0095;
Zw     	=	 0;
Ystar  	=	 0;
Rstar  	=	 0;
PIstar 	=	 0;
q_f =   		 -8;
rf   =  		 0;
f     = 		 0;
Zcp    =		 0;

end;

steady;
resid;
check;

shocks;
//  var etr; stderr SIGMAtr;
//  var ea;  stderr SIGMAa;
//  var er; stderr SIGMAr;
//  var ei; stderr SIGMAi;
//  var epi; stderr SIGMApi;
//  var es; stderr SIGMAs;
//  var emc; stderr SIGMAcp;
    var ew; stderr SIGMAw;

end;
stoch_simul(order=1, irf=30) y, c, mc, i, n, k, rk;
