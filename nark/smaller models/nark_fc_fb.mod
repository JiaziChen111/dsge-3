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

BETA     = 0.999;                // ���η�
DELTA     = 0.025;                // �ں� �����󰢷�
GAMMAtr   = 0.0095;               // ������ �����ݳ� �����
CHIc    = 0.6909;               // �Һ����(consumption habit) ���
KAPPAi    = 3.339;                // �������� ���
CHIn     = 2;                    // parameter on labor function
ALPHAk    = 0.33;                 // �����Լ��� �ں����� ����
THETAh     = 0.807;                // ������ ���� ���� Ȯ�� (1- theta)
THETAf    = 0.546;                // �����籹���� ���� ���� Ȯ�� (1- theta)
THETAw    = 0.541;                // �ӱ� ����Ȯ��
TAYLORpi    = 0.4019;               // ���Ϸ� ��Ģ�� ���÷��̼� ���� ���
TAYLORy    = 0.312;                // ���Ϸ� ��Ģ�� ���갸 ���� ���
TAYLORr    = 0.9407;               // ���Ϸ� ��Ģ�� �ݸ� ���Ӽ� ���
PSIh     = 6.454;                // ������ ��üź�¼�
PSIf    = 5.175;                // ������ ��üź�¼�
PSIw    = 6.091;                // ���� �뵿������ ��üź�¼�   
IOTA    = 0.3162;               // Price indexation ��� 
PHIx   = 0.4880;               // �����Լ� ���Ӽ� ���
XIh  = 0.3945;               // �����Լ��� ȯ�� ź�¼� ���

UIPy     = 0.165;                 // UIP ���ǳ� ���
UIPf     = 0.03;                 // UIP ���ǳ� �ؿܺ�ä�� ���� ź�¼� ���
UIPr    = 0.5;
ETA     = 1.2;                  // �뵿�Լ��� Frisch ���
SP  = 0.975;                // ����� ����Ȯ��
KAPPArp   = 0.05;                 // ����ũ �����̾� �������
ALPHAc   = 0.56723;               // �Һ����� ������ ����
XIc   = 3.955;                // ���ռҺ��� ��ü ź�¼�
ALPHAi   = 0.35088;               // ���ں����� ������ ���� 
XIi   = 1.673;                // ���ں����� ��ü ź�¼�
Gss    = 0.20; 
IT = 1.005;              // ������ǥ���� (1.005^4 = 1.020)   
v       = 1;                    // �����ڱ�(working capital) ���� ����
GAMMA1     = 0.02;                 // capital utilization ���
GAMMA2    = 1.2;                  // capital utilization ���
ALPHAx    = 0.45550;
XIx    = 1.488;
PHIg   = 1;

RHOtr  = 0.8180;                  // ��ݱ��� ���Ӽ� ���
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

SIGMAr   = 0.0014;             // ��� ǥ������
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
                            
a11 =  0.8098;                // Foreign VAR ���  
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
