%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  New Keynesian model with Unit root shock process  %%
%%  2014.4.16                                        %%
%%  Written by Byoungho Bae                 
%%  variable names are modified 2015.2.6 by Yi          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var // detrended variable (cf: Y = Y_trend/Z_t)
    Y, C, N, R, W, Rk, K, dZtr, dP, G 
    // variables required for price stickiness 
    MC, Pstar, X1, X2, SS, ERE
    // variables for open economy
    F, Rf, S, Pxf, Px, Pm, X, PPx, PPm, Rw, Xh, Xf, dPx
    // variables for Financial Accelerator
    Q, I, Re, NW, ut
    // Wage determination
    Wstar, X1w, X1ww,X2w,X2ww,
    // variables for composite goods
    Ch, Cf, P, Pfstar, dPh, dPcf, Pcf, Pif, X1f, X2f, Ih, If, Pi
    Pifstar, dPih, dPif, X1if, X2if, M
    // foreign VAR
    Ystar, Rstar, Pistar
    // shock process 
    Za,  Zi, Zc,  Zpi, Znw, Zpif, Zx, Zg, Zw, Zcp, Zmc,
    // Observation variables
    Y_obs, R_obs,  dP_obs, C_obs, I_obs, NW_obs, dPcf_obs, X_obs, dPi_obs, 
    S_obs, Ystar_obs, Rstar_obs, Pistar_obs, W_obs, G_obs, SP_obs, dPx_obs
    // variables for RBC block (for gap estimation)
    GAPeff, Y_e, Ce, Ne, LAMe, We, R_e,Ke, Inve,  Rke, Qe

    Ystar_or;

varexo er, 
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
parameters  BETA, DELTA, CHIn, ALPHAk,  CHIc, PSIh, PSIf, THETAh, THETAf, TAYLORpi, TAYLORy, 
            TAYLORr, UIPy, UIPf, PHIx, XIh, KAPPAi, GAMMAtr, SP, KAPPArp, 
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
///////////////////////////// �Һ������ //////////////////////////////////
exp(dZtr(+1))* exp(Zc)*1/(exp(C)-CHIc*exp(C(-1))/exp(dZtr)) = 
    BETA*exp(Zc(+1))*1/(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1))) *exp(R)/exp(dP(+1));   


/////////////  �ӱݰ����� ���� ������ /////////////////  
exp(X1w) = exp(Wstar)^(-PSIw*(ETA+1))*exp(X1ww);
exp(X1ww) = exp(W)^(PSIw*(ETA+1))*exp(N)^(ETA+1) 
            + BETA*THETAw*exp(dZtr(+1))^(PSIw*(ETA+1))*exp(X1ww(+1));
exp(X2w) = exp(Wstar)^(-PSIw)*exp(X2ww);
exp(X2ww) = 1/(exp(C)-CHIc*exp(C(-1))/exp(dZtr))*exp(W)^PSIw*exp(N) 
            + BETA*THETAw*exp(dZtr(+1))^(PSIw-1)*exp(X2ww(+1));
exp(Wstar)= PSIw/(PSIw-1)*exp(Zw)*exp(X1w) / exp(X2w);
exp(W)   = (THETAw*exp(W(-1)/exp(dZtr))^(1-PSIw) 
            + (1-THETAw)*exp(Wstar)^(1-PSIw))^(1/(1-PSIw));                                                     

/////////////////////////// �Һ����� ���� ���ǽ� /////////////////////////
exp(Ch) = ALPHAc*exp(P)^XIc*exp(C);                        
exp(Cf) = ((1-ALPHAc)*(exp(P)/exp(Pcf))^XIc*exp(C));                
exp(P) = (ALPHAc + (1-ALPHAc)*(exp(Pcf))^(1-XIc))^(1/(1-XIc));   
exp(dP)*exp(Zpi) = (exp(P)/exp(P(-1))*exp(dPh));                

///////////////////////////// ���� ������ ���ǽ� (���� �� ���� �Һ���) //////
exp(X1) = exp(Y)*exp(MC) +
          THETAh*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr(-1)))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr)))
            /exp(dP(+1))*exp(dPh(+1))^(Zmc+1)*exp(dPh)^((IOTA)*(-Zmc))
            *IT^((-Zmc)*(1-IOTA))*exp(X1(+1));
exp(X2) = exp(Y) +
          THETAh*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr(-1)))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr)))
            /exp(dP(+1))*exp(dPh(+1))^(Zmc)*exp(dPh)^((IOTA)*(1-Zmc))*IT^((1-Zmc)
            *(1-IOTA))*exp(X2(+1));
exp(Pstar) = Zmc/(Zmc-1)* exp(X1)/exp(X2);
1 = exp(dPh(-1))^(IOTA*(1-Zmc))*IT^((1-IOTA)*(1-Zmc))*THETAh*exp(dPh)^(Zmc-1) +
    (1-THETAh)*exp(Pstar)^(1-Zmc);
exp(SS) = (exp(dPh(-1))^IOTA*IT^(1-IOTA))^(-PSIh)*THETAh*exp(dPh)^(PSIh)*exp(SS(-1)) +
          (1-THETAh)*exp(Pstar)^(-PSIh);

exp(X1f) = exp(Cf)*exp(ERE)*exp(P)*(1+v*(exp(R)-1))*exp(Pistar) + 
           THETAf*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr(-1)))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr)))
            /exp(dP(+1))*(exp(dPcf(+1)))^(PSIf)*(exp(dPcf))^((IOTA)*(-PSIf))
            *IT^(-PSIf*(1-IOTA))*exp(dPh(+1))*exp(X1f(+1));
exp(X2f) = exp(Cf) + 
           THETAf*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr(-1)))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr)))
            /exp(dP(+1))*(exp(dPcf(+1)))^(PSIf)*(exp(dPcf))^((IOTA)*(1-PSIf))
            *IT^((1-PSIf)*(1-IOTA))*exp(X2f(+1));
exp(Pfstar) = PSIf/(PSIf-1)* exp(X1f)/ exp(X2f);
1 = exp(dPcf(-1))^(IOTA*(1-PSIf))*IT^((1-IOTA)*(1-PSIf))* THETAf*exp(dPcf)^(PSIf-1) + 
    (1-THETAf)*(exp(Pfstar)/exp(Pcf))^(1-PSIf);
exp(dPcf)*exp(Zpif) = exp(Pcf)/exp(Pcf(-1))*exp(dPh);
 
/////////////////////////// ���ں����� ���� ���ǽ� ////////////////////////
exp(Ih) = ALPHAi*exp(Pi)^XIi*exp(I);
exp(If) = (1-ALPHAi)*(exp(Pi)/exp(Pif))^XIi*exp(I);
exp(Pi) = (ALPHAi + (1-ALPHAi)*(exp(Pif))^(1-XIi))^(1/(1-XIi));
exp(dPih) = exp(Pi)/exp(Pi(-1))*exp(dPh);

exp(X1if) = exp(If)*exp(ERE)*exp(P)*(1+v*(exp(R)-1))*exp(Pistar) + 
               THETAf*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr(-1)))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr)))
                /exp(dP(+1))*(exp(dPif(+1)))^(PSIf)*(exp(dPif))^((IOTA)*(-PSIf))
                *IT^(-PSIf*(1-IOTA))*exp(dPh(+1))*exp(X1if(+1));
exp(X2if) = exp(If) + 
               THETAf*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr(-1)))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr)))
                /exp(dP(+1))*(exp(dPif(+1)))^(PSIf)*(exp(dPif))^((IOTA)*(1-PSIf))
                *IT^((1-PSIf)*(1-IOTA))*exp(X2if(+1));
exp(Pifstar) = PSIf/(PSIf-1)* exp(X1if)/ exp(X2if);
1 = exp(dPif(-1))^(IOTA*(1-PSIf))*IT^((1-IOTA)*(1-PSIf))*THETAf*(exp(dPif))^(PSIf-1) + 
    (1-THETAf)*(exp(Pifstar)/exp(Pif))^(1-PSIf);
exp(dPif)*exp(Zpif) = exp(Pif)/exp(Pif(-1))*exp(dPh);

///////////////// �������ӱ� (Financial Accelerator) ///////////////////////
exp(Re) = (exp(Q)*(1-DELTA)+ exp(Rk)*exp(ut)-1/exp(Zi)*(GAMMA1*(exp(ut)-1)+GAMMA2/2*(exp(ut)-1)^2))
            *exp(dPh)/exp(Q(-1));
exp(NW) = exp(Znw)*SP*(exp(Re)*exp(Q(-1))*exp(K(-1))/exp(dZtr)-
                         (exp(NW(-1))/(exp(K(-1))*exp(Q(-1))))^(-KAPPArp)*exp(R(-1))/exp(dP)
                           *(exp(Q(-1))*exp(K(-1))/exp(dZtr)-exp(NW(-1))/exp(dZtr)));
exp(Re(+1)) = (exp(NW)/(exp(K)*exp(Q)))^(-KAPPArp)*exp(R)/exp(dP(+1));

%% ���� ���½� 
exp(K)-(1-DELTA)*exp(K(-1))/exp(dZtr) = exp(Zi)*exp(I)
                                      *(1-KAPPAi/2*(exp(I)*exp(dZtr)/exp(I(-1))-exp(GAMMAtr))^2);

%% �ں� �������� ���� �ϰ����ǽ� Rk = A(ut)'
exp(Rk) = (GAMMA1+GAMMA2*(exp(ut)-1))/exp(Zi);

%% ���ڿ� ���� �ϰ����ǽ� (�ں��� ���� ������)                                                                                              
1 = exp(Q)*exp(Zi)*(1 - 
                      KAPPAi/2*(exp(I)*exp(dZtr)/exp(I(-1)) -exp(GAMMAtr))^2 -
                      KAPPAi*(exp(I)*exp(dZtr)/exp(I(-1))-exp(GAMMAtr))
                        *exp(I)*exp(dZtr)/exp(I(-1))) + 
    BETA*exp(Q(+1))*exp(Zi(+1))
      *(1/exp(dZtr(+1))*(exp(C)-CHIc*exp(C(-1))/exp(dZtr))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1))))
        *KAPPAi*(exp(I(+1))*exp(dZtr(+1))/exp(I)-exp(GAMMAtr))
          *(exp(I(+1))*exp(dZtr(+1))/exp(I))^2;


////////////////////////////////// ����ι� ///////////////////////////////
// �ں� �뵿�� �������Կ�� ����
exp(W)*(1+v*(exp(R)-1))/exp(Rk) = (1-ALPHAk)/ALPHAk*exp(K(-1))*exp(ut)/exp(N)*1/exp(dZtr);

// ��ũ�� ���ǽ�
exp(MC) = 1/exp(Za)*(exp(W)*(1+v*(exp(R)-1)))^(1-ALPHAk)*exp(Rk)^ALPHAk
            *(1-ALPHAk)^(ALPHAk-1)*ALPHAk^(-ALPHAk);

//exp(W)*(1+v*(exp(R)-1)) = exp(Za)*exp(MC)*(1-ALPHAk)*(exp(K(-1))*exp(ut))^(ALPHAk)
//                            *exp(N)^(-ALPHAk)*exp(dZtr)^(-ALPHAk)/exp(SS);
//exp(Rk) = exp(Za)*exp(MC)*ALPHAk*(exp(K(-1))*exp(ut))^(ALPHAk-1)
//            *exp(N)^(1-ALPHAk)*exp(dZtr)^(1-ALPHAk)/exp(SS);

exp(Y) = exp(Za)*(exp(K(-1))*exp(ut))^ALPHAk*exp(N)^(1-ALPHAk)*exp(dZtr)^(-ALPHAk)/exp(SS);
exp(Y) = exp(Ch) + exp(Ih) + exp(Xh) + exp(G) +
         KAPPAi/2*(exp(I)*exp(dZtr)/exp(I(-1))-exp(GAMMAtr))^2*exp(I) + 
         1/exp(Zi)*(GAMMA1*(exp(ut)-1)+GAMMA2*(exp(ut)-1)^2)*exp(ut)*exp(K(-1));
exp(G) = Gss*exp((Y))*exp(Zg)*(exp(Y)/exp(steady_state(Y)))^(-PHIg);


/////////////////////////////////// ������� ���μ��� //////////////////////
// ������ ������ ln(A_t) = ln(A_t-1) + GAMMAtr*exp(Zg) + ea //
dZtr         = RHOtr*dZtr(-1) + (1-RHOtr)*GAMMAtr + etr;                                                
Za         = RHOa    *Za(-1)      + ea;
Zi       = RHOi  *Zi(-1)    + ei;
Zc         = RHOc    *Zc(-1)      + ec;
Zpi       = RHOpi  *Zpi(-1)    + epi;
Znw         = RHOnw    *Znw(-1)      + enw;
Zpif      = RHOpif *Zpif(-1)   + epif;
Zx         = RHOx    *Zx(-1)      + ex;
Zg         = RHOg    *Zg(-1)      + eg;
Zw         = RHOw    *Zw(-1)      + ew;
Zcp        = RHOcp   *Zcp(-1)     + ecp;
Zmc        = (1-RHOmc)*PSIh + RHOmc   *Zmc(-1)   + emc;


/////////////////////////// ��ȭ��Ģ (Taylor rule) ////////////////////////
log(exp(R)) = TAYLORr*log(exp(R(-1))) + 
              (1-TAYLORr)*(log(steady_state(exp(R)))+(1+TAYLORpi)*log(exp(dP)) +
                        TAYLORy*log(exp(Y)/steady_state(exp(Y)))) + 
              er;


////////////////////////////// Foreign VAR/////////////////////////////////
Ystar   = a11*Ystar(-1)  + a12*Rstar(-1) + a13*Pistar(-1) + eystar;
Rstar   = a21*Ystar(-1)  + a22*Rstar(-1) + a23*Pistar(-1) + erstar;
Pistar  = a31*Ystar(-1)  + a32*Rstar(-1) + a33*Pistar(-1) + epistar;

 
/////////////////////////// Export function & BOP  ////////////////////////
// �����Լ�
exp(X) = exp(Zx)*(exp(X(-1))/exp(dZtr))^PHIx
                  *(((1+v*(exp(R)-1))*exp(Px)/exp(S))^(-XIh)*(exp(Ystar)))^(1-PHIx);
exp(Xh) = ALPHAx*exp(Px)^XIx*exp(X);
exp(Xf) = (1-ALPHAx)*(exp(Px)/(exp(ERE)*exp(P)))^XIx*exp(X);
exp(Px) = (ALPHAx+(1-ALPHAx)*(exp(ERE)*exp(P))^(1-XIx))^(1/(1-XIx));
exp(dPx)*exp(epx) = exp(Px)/exp(Px(-1))*exp(dPh);

exp(ERE)*exp(Pxf)*exp(P)*exp(X) + 
exp(S)*exp(Rf(-1))*exp(F(-1))*exp(P)/exp(dZtr)/exp(dP) 
= exp(S)*exp(F)*exp(P) + exp(ERE)*exp(P)*exp(Cf)*(1+v*(exp(R)-1)) + 
  exp(ERE)*exp(P)*exp(If)*(1+v*(exp(R)-1)) + 
  exp(ERE)*exp(P)*exp(Xf)*(1+v*(exp(R)-1));

//////////////////////////////  UIP condition /////////////////////////////
// UIP condition
exp(R) = exp(Rf)*exp(S(+1))/(exp(S));
exp(Rf) = exp(-UIPf*(exp(S)*exp(F)- (exp(steady_state(Y))*UIPy)) - 
              UIPr*(exp(Rf)-exp(steady_state(Rf))- (exp(R)-exp(steady_state(R)))))
            *exp(Zcp)*exp(Rw)*(exp(Rstar));

exp(Rw) = (1+GAMMAtr)/BETA;
exp(Pxf)  = exp(PPx);//*exp(S));                                         
exp(Pm)  = exp(PPm)*(exp(Pistar))*exp(S);            // ������ ����
exp(ERE) = exp(PPm)*(exp(Pistar))*exp(S)/exp(P);    // ����ȯ��
exp(PPx) = 1;
exp(PPm) = 1;
exp(M) = exp(Cf) + exp(If) + exp(Xf);
  

/////////// Efficient Gap ������ ���� RBC block /////////////////////////
// RBC (1)
exp(LAMe) = 1/(exp(Ce));

// RBC (2)
CHIn/(1-exp(Ne)) = 1/exp(Ce)*exp(We);    //CHIn*exp(Ne)^ETA = exp(LAMe)*exp(We);

// RBC (3)
exp(dZtr(+1))*exp(Qe)*1/exp(Ce) = BETA*1/exp(Ce(+1))*(exp(Rke(+1))+exp(Qe(+1))*(1-DELTA));

// RBC (4)
exp(dZtr(+1))*1/exp(Ce) = BETA*(1/exp(Ce(+1))*exp(R_e));

// RBC (5)
exp(We) = (1-ALPHAk)*exp(Za)*exp(Ke(-1))^(ALPHAk)*exp(Ne)^(-ALPHAk) *exp(dZtr)^(1-ALPHAk);

// RBC (6)
exp(Rke) = ALPHAk*exp(Za)*exp(Ke(-1))^(ALPHAk-1)*exp(Ne)^(1-ALPHAk)*exp(dZtr)^(1-ALPHAk);

// RBC (7)
exp(Y_e) = exp(Za)*exp(Ke(-1))^ALPHAk*exp(Ne)^(1-ALPHAk)*exp(dZtr)^(-ALPHAk);

// RBC (8)
exp(Y_e) = exp(Ce) + exp(Ke) - (1-DELTA)*exp(Ke(-1))*1/exp(dZtr) + 
           KAPPAi/2*(exp(Inve)/exp(Ke)-DELTA)^2*exp(Ke);

// RBC (9)
exp(Qe) = (1+KAPPAi*(exp(Inve)/ exp(Ke) - DELTA));

// RBC (10)
exp(Inve)*exp(Zi) = exp(Ke)- (1- DELTA)*exp(Ke(-1))*1/exp(dZtr);

GAPeff = (exp(Y)-exp(steady_state(Y)))/exp(steady_state(Y)) - 
         (exp(Y_e)-exp(steady_state(Y_e)))/exp(steady_state(Y_e));

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


///////////////////////// ����������(Measurement Equations) ////////////////
Y_obs       =   100*(Y-Y(-1) + dZtr) - 100*GAMMAtr;
R_obs       =   400*(R - steady_state(R));                      
C_obs       =   100*(C-C(-1) + dZtr) - 100*GAMMAtr;
dP_obs      =   100*(dP - steady_state(dP)); 
I_obs     =   100*(I-I(-1) + dZtr) - 100*GAMMAtr;
G_obs       =   100*(G-G(-1) + dZtr) - 100*GAMMAtr;
NW_obs      =   100*(NW-NW(-1));  
dPcf_obs     =   100*(dPcf) + 0.94;                 
dPx_obs     =   100*(dPx) - 0.258;                         
X_obs       =   100*(X-X(-1) + dZtr) - 100*GAMMAtr;
S_obs       =   100*(ERE-ERE(-1)); 
dPi_obs   =   400*(dPih-steady_state(dPih));
SP_obs      =   400*(Re-R - (steady_state(Re) - steady_state(R)))                           ;  
Ystar_obs   =   100*(Ystar-Ystar(-1));
Rstar_obs   =   400*(Rstar); // steady state value or Rstar_obs  = 2.88, based on data
Pistar_obs  =   400*(Pistar); // steady state value or Rstar_obs  = 2.17. based on data
W_obs       =   100*(W-W(-1) + dZtr) - 100*GAMMAtr;
/////////////////////////////////////////////////////////////////////

Ystar_or = (Y-Y(-4))*100;


end;

observation_trends;
Y_obs (GAMMAtr);
C_obs (GAMMAtr);
I_obs (GAMMAtr);
X_obs (GAMMAtr);
G_obs (GAMMAtr);
W_obs (GAMMAtr);

end;


initval;
Y           = 	 0.572716;
C           =    0.0854136;
N           =	-0.184133;
R          	=	 log((1+GAMMAtr)/BETA);//0.0125045;
W          	=	 0.161546;
Rk         	=	-2.82762;
K          	=	 2.11885;
dZtr         	=	 GAMMAtr;
dP         	=	 0;
G          	=	-1.03672;
MC         	=	-0.182322;
Pstar          	=	 0;
X1         	=	 1.74009;
X2         	=	 1.92241;
SS         	=	 0;
ERE        	=	-0.286734;
F          	=	-0.227666;
Rf         	=	 0.0125045;
S          	=	-0.286735;
Pxf         	=	 1.7212e-06;
Pm         	=	-0.286734;
X          	=	-0.631171;
PPx        	=	 4.00498e-07;
PPm        	=	 4.00498e-07;
Rw        	=	 0.0124597;
Q          	=	 2.23411e-11;
I        	=	-1.25595;
Re         	=	 0.0348178;
NW         	=	 1.89572;
Ch         	=	-0.24304;
Cf         	=	-1.16408;
P         	=	-0.0351037;
Pfstar         	=	-0.155084;
dPh        	=	  0;
dPcf        	=	-0.0990233;
X1f        	=	-0.97575;
X2f        	=	-0.638344;
Ih       	=	-2.65812;
If       	=	-1.5351;
Pi       	=	-0.06606;
Pifstar     	=	-0.0919085;
dPih      	=	 0;
dPif    	=	-0.0432459;
X1if    	=	-0.0223043;
X2if    	=	-0.0223043;
M        	=	-0.639335;
Zmc        	=	 6;
GAPeff     	=	 0;
Y_e        	=	-0.0751586;
Ce         	=	-0.423288;
Ne         	=	-1.12738;
LAMe       	=	 0.423288;
We         	=	 0.661244;
R_e        	=	 0.0125075;
Ke         	=	 2.07067;
Inve       	=	-1.30431;
Rke        	=	-3.24499;
Qe         	=	 0.036211;
Wstar          	=	 0.161546;
X1w        	=	 1.19985;
ut         	=	 0.0632151;
Xh          =    log(0.2);
Xf          =    log(0.2);
Px         =    log(1);
end;

steady;
resid;
check;

shocks;
//  var etr; stderr SIGMAtr;
//  var ea; stderr SIGMAa;
 // var er; stderr SIGMAr;
//  var ei; stderr SIGMAi;
//  var epi; stderr SIGMApi;
//  var es; stderr SIGMAs;
  var ecp; stderr SIGMAcp;
//  var ew; stderr SIGMAw;

end;
stoch_simul(order=1, irf=30) Y, C, I, dP, R, X,S, F;

//write_latex_dynamic_model;

/////////////////////////////////////////
//      Bayesian Estimation            //
//    (2001. 1/4~ 2012. 3/4)           //
/////////////////////////////////////////  

varobs Y_obs, R_obs, dP_obs, C_obs, I_obs, NW_obs, X_obs,S_obs, W_obs, Ystar_obs, 
Rstar_obs, Pistar_obs, G_obs, dPx_obs;
//, dPcf_obs;//,  G_obs;//;// dPcf_obs;//, W_obs ;//dPcf_obs,  dPi_obs ;
//,  Ystar_obs, Rstar_obs, Pistar_obs G_obs, W_obs, dPd_obs;//, N_obs;//, SP_obs;

estimated_params;
    ETA,        gamma_pdf,      1.4,     0.2;
    PSIh,        gamma_pdf,      6,       1.0;
    PSIf,       gamma_pdf,      6,       1.0;
    PSIw,       gamma_pdf,      6,       1.0;
    THETAh,        beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05 
    THETAf,       beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05
    THETAw,       beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05
    CHIc,       normal_pdf,     0.5,     0.1;
    KAPPAi,       normal_pdf,     4,       1.0;
    KAPPArp,      beta_pdf,       0.3,     0.05;    
    
    TAYLORr,       beta_pdf,       0.8,     0.05;
    TAYLORy,       beta_pdf,       0.2,     0.05;
    TAYLORpi,       beta_pdf,       0.5,     0.1;   
    IOTA,       beta_pdf,       0.5,     0.1;
    PHIx,      beta_pdf,       0.5,     0.1;
    XIh,     gamma_pdf,      0.5,     0.1;
    ALPHAc,      beta_pdf,       0.6,     0.1;
    ALPHAi,      beta_pdf,       0.6,     0.1;     
    XIc,      gamma_pdf,      2,       0.5;
    XIi,      gamma_pdf,      2,       0.5;

    ALPHAx,       beta_pdf,       0.5,     0.1;
    XIx,       gamma_pdf,      2.5,     0.5;
    GAMMA2,       beta_pdf,       0.5,     0.1;     
    UIPr,       beta_pdf,       0.5,     0.1;     
    PHIg,      gamma_pdf,      0.5,     0.1;
    
    RHOa,       beta_pdf,       0.8,    0.1;
    RHOc,       beta_pdf,       0.8,    0.1;
    RHOtr,     beta_pdf,       0.8,    0.1;
    RHOpi,     beta_pdf,       0.8,    0.1;
    RHOi,     beta_pdf,       0.8,    0.1;
    RHOnw,       beta_pdf,       0.8,    0.1;
    RHOpif,    beta_pdf,       0.8,    0.1;
    RHOx,       beta_pdf,       0.8,    0.1;
    RHOpii,    beta_pdf,       0.8,    0.1;
    RHOg,       beta_pdf,       0.8,    0.1;
    RHOw,       beta_pdf,       0.8,    0.1;
    RHOcp,      beta_pdf,       0.8,    0.1;

    a11,        normal_pdf,     0.9,    0.1;
    a12,        normal_pdf,     0.0,    0.1;
    a13,        normal_pdf,     0.0,    0.1;
    a21,        normal_pdf,     0.0,    0.1;
    a22,        normal_pdf,     0.9,    0.1;
    a23,        normal_pdf,     0.0,    0.1;
    a31,        normal_pdf,     0.0,    0.1;
    a32,        normal_pdf,     0.0,    0.1;
    a33,        normal_pdf,     0.6,    0.1;

    stderr er,      inv_gamma_pdf, 0.01,  2;
    stderr ea,      inv_gamma_pdf, 0.01,  2;
    stderr ec,      inv_gamma_pdf, 0.01,  2;
    stderr etr,  inv_gamma_pdf, 0.01,  2;
    stderr epi,    inv_gamma_pdf, 0.01,  2;
    stderr ei,    inv_gamma_pdf, 0.01,  2;
    stderr enw,      inv_gamma_pdf, 0.01,  2;
    stderr epif,   inv_gamma_pdf, 0.01,  2;
    stderr ex,      inv_gamma_pdf, 0.01,  2;
    stderr eg,      inv_gamma_pdf, 0.01,  2;
    stderr emc,     inv_gamma_pdf, 0.01,  2;
    stderr eystar,  inv_gamma_pdf, 0.01,  2;
    stderr erstar,  inv_gamma_pdf, 0.01,  2;
    stderr epistar, inv_gamma_pdf, 0.01,  2;
    stderr ew,      inv_gamma_pdf, 0.01,  2;
    stderr ecp,     inv_gamma_pdf, 0.01,  2;
    stderr epx,     inv_gamma_pdf, 0.01,  2;



end;

/*
//estimation(
datafile     = data_201414, //data_201344_sa_soe_dmean.xls, //,
first_obs    = 1, 
mh_replic    = 50000,
//mode_file    = BOKDSGE2014_0623_TR_UIP_W_markup_mode,
mh_nblocks   = 1,
//forecast     = 12, 
//nobs       = [5:50],
mh_drop      = 0.4, 
mh_jscale    = 0.20, 
mode_compute = 4, 
plot_priors  = 0,
smoother, 
//endogenous_prior,
//,mode_check
filtered_vars,
nograph) Y_obs, R_obs,  dP_obs, C_obs, I_obs, NW_obs, dPcf_obs, X_obs, dPi_obs, 
S_obs, Ystar_obs, Rstar_obs, Pistar_obs dZtr G_obs W_obs Za,  Zi, Zc,  
Zpi, Znw, Zpif, Zx, Zg, Zw, Zcp;

shock_decomposition   Y_obs, C_obs, R_obs, dP_obs, I_obs, GAPeff;
*/   

            
     