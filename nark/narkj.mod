%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  New Keynesian model with Unit root shock process  %%
%%  2014.4.16                                        %%
%%  Written by Byoungho Bae                 
%%  variable names are modified 2015.2.6 by Yi          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var // detrended variable (cf: Y = Y_trend/Z_t)
    Y, C, N, R, W, Rk, K, dZtr, dP, G 
    // variables required for price stickiness 
    MC, Phstar, X1, X2, SS, Q
    // variables for open economy
    F, Rf, S, Pxf, Px, Pm, X, PPx, PPm, Rw, Xh, Xf, dPx
    // variables for Financial Accelerator
    Pk, I, Re, NW, ut
    // Wage determination
    Wstar, X1w, X1ww,X2w,X2ww,
    // variables for composite goods
    Ch, Cf, P, Pcfstar, dPh, dPcf, Pcf, Pif, X1f, X2f, Ih, If, Pi
    Pifstar, dPi, dPif, X1if, X2if, M
    // foreign VAR
    Ystar, Rstar, PIstar
    // shock process 
    Za,  Zi, Zc,  Zpi, Znw, Zpif, Zx, Zg, Zw, Zcp, Zmc, Zq
    // Observation variables
    Y_obs, R_obs,  dP_obs, C_obs, Inv_obs, NN_obs, dPcf_obs, X_obs, dPi_obs, 
    S_obs, Ystar_obs, Rstar_obs, Pistar_obs, G_obs, SP_obs, dPx_obs,
    //W_obs
    // variables for RBC block (for gap estimation)
    // GAPeff, Y_e, Ce, Ne, LAMe, We, R_e,Ke, Inve,  Rke, Qe, 
    dY

    ///* fb F-Blocks
    c, r, w, n, ch, cf, pcf, ph, mc, mcf,
    q, ih, if, pi, i, pif, re, pk, rk, u,
    nw, k, y, xh, g, x, px, pxf, xf, rf,
    f, m,
    Ygap, Rgap
    //fe

    // Flexible-Closed-No FA block: 
    ct, rt, wt, kt, nt, mct, yt, rkt, it, 
    ret, nwt, pkt, utt,
    Ygapt, Rgapt

;

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
       epx,
       eq;
       

////////////// define parameters //////////////
parameters  BETA, DELTA, CHIn, ALPHAk,  CHIc, PSIh, PSIf, THETAh, THETAf, TAYLORpi, 
            TAYLORy, TAYLORr, UIPy, UIPf, PHIx, XIh, KAPPAi, GAMMAtr, SP, KAPPArp, 
            ALPHAc, XIc, ALPHAi, XIi, ETA, IOTA, IT, THETAw, 
            PSIw, v, GAMMA1, GAMMA2, UIPr, ALPHAx, XIx, RHOtr, RHOa, RHOi, 
            RHOc,  RHOpi, RHOnw, RHOpif, RHOx, RHOpii, Gss, RHOg,
            RHOw, RHOcp, RHOmc, RHOq, SIGMAr, SIGMAa, SIGMAi, SIGMAc,  SIGMApi,
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
RHOmc   = 0.5;          // persistence of shock on marginal cost
RHOq    = 0.3;          // persistence of shock on the relation between r, rf and q

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
///////////////////////////// Euler //////////////////////////////////
///////////////////////////////////////////////////////////////////////////
exp(dZtr(+1))*exp(Zc)/(exp(C)-CHIc*exp(C(-1))/exp(dZtr))
= BETA*exp(Zc(+1))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1)))*exp(R)/exp(dP(+1));   


////////////////////////////  W ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
exp(X1w) = exp(Wstar)^(-PSIw*(ETA+1))*exp(X1ww);
exp(X1ww) = exp(W)^(PSIw*(ETA+1))*exp(N)^(ETA+1) 
            + BETA*THETAw*exp(dZtr(+1))^(PSIw*(ETA+1))*exp(X1ww(+1));
exp(X2w) = exp(Wstar)^(-PSIw)*exp(X2ww);
exp(X2ww) = 1/(exp(C)-CHIc*exp(C(-1))/exp(dZtr))*exp(W)^PSIw*exp(N) 
            + BETA*THETAw*exp(dZtr(+1))^(PSIw-1)*exp(X2ww(+1));
exp(Wstar)= PSIw/(PSIw-1)*exp(Zw)*exp(X1w)/exp(X2w);
exp(W) = (THETAw*exp(W(-1)/exp(dZtr))^(1-PSIw)
         + (1-THETAw)*exp(Wstar)^(1-PSIw))^(1/(1-PSIw));


/////////////////////// C and P /////////////////////////////////
/////////////////////////////////////////////////////////////////
exp(Ch) = ALPHAc*exp(P)^XIc*exp(C);
exp(Cf) = ((1-ALPHAc)*(exp(P)/exp(Pcf))^XIc*exp(C));
exp(P) = (ALPHAc+(1-ALPHAc)*(exp(Pcf))^(1-XIc))^(1/(1-XIc));
exp(dP)*exp(Zpi) = (exp(P)/exp(P(-1))*exp(dPh));


////////////////// Ph ////////////////////////////
///////////////////////////////////////////////////////////////////////////
exp(X1) = exp(Y)*exp(MC) 
          + THETAh*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr))
                         /(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1)))) //?? dZtr
            /exp(dP(+1))*exp(dPh(+1))^(PSIh+1)*exp(dPh)^((IOTA)*(-PSIh))
            *IT^((-PSIh)*(1-IOTA))*exp(X1(+1));
exp(X2) = exp(Y)
          + THETAh*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr))
                         /(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1)))) //?? dZtr
            /exp(dP(+1))*exp(dPh(+1))^PSIh*exp(dPh)^((IOTA)*(1-PSIh))
            *IT^((1-PSIh)*(1-IOTA))*exp(X2(+1));
exp(Phstar) = PSIh/(PSIh-1)*exp(Zmc)*exp(X1)/exp(X2);
1 = exp(dPh(-1))^(IOTA*(1-PSIh))*IT^((1-IOTA)*(1-PSIh))*THETAh*exp(dPh)^(PSIh-1)
    + (1-THETAh)*exp(Phstar)^(1-PSIh);
exp(SS) = (exp(dPh(-1))^IOTA*IT^(1-IOTA))^(-PSIh)*THETAh*exp(dPh)^(PSIh)*exp(SS(-1))
          + (1-THETAh)*exp(Phstar)^(-PSIh);


////////////////////// Pcf ////////////////////////////
///////////////////////////////////////////////////////////////////////////
exp(X1f) = exp(Cf)*exp(Q)*exp(P)*(1+v*(exp(R)-1))*exp(PIstar)
           + THETAf*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr))
                          /(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1)))) //?? dZtr
             /exp(dP(+1))*(exp(dPcf(+1)))^(PSIf)*(exp(dPcf))^((IOTA)*(-PSIf))
             *IT^(-PSIf*(1-IOTA))*exp(dPh(+1))*exp(X1f(+1));
exp(X2f) = exp(Cf)
           + THETAf*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr))
                          /(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1)))) //?? dZtr
             /exp(dP(+1))*(exp(dPcf(+1)))^(PSIf)*(exp(dPcf))^((IOTA)*(1-PSIf))
             *IT^((1-PSIf)*(1-IOTA))*exp(X2f(+1));
exp(Pcfstar) = PSIf/(PSIf-1)*exp(X1f)/exp(X2f);
1 = exp(dPcf(-1))^(IOTA*(1-PSIf))*IT^((1-IOTA)*(1-PSIf))*THETAf*exp(dPcf)^(PSIf-1)
    + (1-THETAf)*(exp(Pcfstar)/exp(Pcf))^(1-PSIf);
exp(dPcf)*exp(Zpif) = exp(Pcf)/exp(Pcf(-1))*exp(dPh);
 
/////////////////////////// I and Pi ////////////////////////
///////////////////////////////////////////////////////////////////////////
exp(Ih) = ALPHAi*exp(Pi)^XIi*exp(I);
exp(If) = (1-ALPHAi)*(exp(Pi)/exp(Pif))^XIi*exp(I);
exp(Pi) = (ALPHAi+(1-ALPHAi)*(exp(Pif))^(1-XIi))^(1/(1-XIi));
exp(dPi) = exp(Pi)/exp(Pi(-1))*exp(dPh);

exp(X1if) = exp(If)*exp(Q)*exp(P)*(1+v*(exp(R)-1))*exp(PIstar)
            + THETAf*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr))
                           /(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1)))) ///// dZtr
              /exp(dP(+1))*(exp(dPif(+1)))^(PSIf)*(exp(dPif))^((IOTA)*(-PSIf))
              *IT^(-PSIf*(1-IOTA))*exp(dPh(+1))*exp(X1if(+1));
exp(X2if) = exp(If)
            + THETAf*BETA*((exp(C)-CHIc*exp(C(-1))/exp(dZtr))
                           /(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1)))) ///// dZtr
              /exp(dP(+1))*(exp(dPif(+1)))^(PSIf)*(exp(dPif))^((IOTA)*(1-PSIf))
              *IT^((1-PSIf)*(1-IOTA))*exp(X2if(+1));
exp(Pifstar) = PSIf/(PSIf-1)*exp(X1if)/exp(X2if);
1 = exp(dPif(-1))^(IOTA*(1-PSIf))*IT^((1-IOTA)*(1-PSIf))*THETAf*(exp(dPif))^(PSIf-1)
    + (1-THETAf)*(exp(Pifstar)/exp(Pif))^(1-PSIf);
exp(dPif)*exp(Zpif) = exp(Pif)/exp(Pif(-1))*exp(dPh);


/////////////////////////// Export function & BOP  ////////////////////////
///////////////////////////////////////////////////////////////////////////
//?? no foreign price?
exp(X) = exp(Zx)*(exp(X(-1))/exp(dZtr))^PHIx
         *(((1+v*(exp(R)-1))*exp(Px)/exp(S))^(-XIh)*(exp(Ystar)))^(1-PHIx); 
//exp(X) = exp(Zx)*(exp(X(-1))/exp(dZtr))^PHIx
//         *(((1+v*(exp(R)-1))*exp(Px)/exp(S)/exp(PIstar))^(-XIh)*(exp(Ystar)))^(1-PHIx);
exp(Xh) = ALPHAx*exp(Px)^XIx*exp(X);
exp(Xf) = (1-ALPHAx)*(exp(Px)/(exp(Q)*exp(P)))^XIx*exp(X);
exp(Px) = (ALPHAx+(1-ALPHAx)*(exp(Q)*exp(P))^(1-XIx))^(1/(1-XIx));
exp(dPx)*exp(epx) = exp(Px)/exp(Px(-1))*exp(dPh);

//exp(Px)*exp(X) //?? exp(Q)*exp(Pxf)*exp(P)*exp(X) //this modification chages SS values
exp(Q)*exp(Pxf)*exp(P)*exp(X)
+ exp(S)*exp(Rf(-1))*exp(F(-1))*exp(P)/exp(dZtr)/exp(dP) 
= exp(S)*exp(F)*exp(P)
  + exp(Q)*exp(P)*(1+v*(exp(R)-1))*(exp(Cf) + exp(If) + exp(Xf));


///////////////// Financial Accelerator //////////////////////
//////////////////////////////////////////////////////////////
exp(Re) = (exp(Pk)*(1-DELTA) + exp(Rk)*exp(ut)
           - 1/exp(Zi)*(GAMMA1*(exp(ut)-1)+GAMMA2/2*(exp(ut)-1)^2))
          *exp(dPh)/exp(Pk(-1));
exp(NW) = exp(Znw)*SP
          *(exp(Re)*exp(Pk(-1))*exp(K(-1))/exp(dZtr)
            - (exp(NW(-1))/(exp(K(-1))*exp(Pk(-1))))^(-KAPPArp)*exp(R(-1))/exp(dP)
              *(exp(Pk(-1))*exp(K(-1))/exp(dZtr)-exp(NW(-1))/exp(dZtr)));
exp(Re(+1)) = (exp(NW)/(exp(K)*exp(Pk)))^(-KAPPArp)*exp(R)/exp(dP(+1));

// Law of Motion for Capital
exp(K)-(1-DELTA)*exp(K(-1))/exp(dZtr)
= exp(Zi)*exp(I)*(1-KAPPAi/2*(exp(I)*exp(dZtr)/exp(I(-1))-exp(GAMMAtr))^2);

// Rk = A(ut)'
exp(Rk) = (GAMMA1+GAMMA2*(exp(ut)-1))/exp(Zi);

// FOC for Investment, Pk
1 = exp(Pk)*exp(Zi)*(1 - KAPPAi/2*(exp(I)*exp(dZtr)/exp(I(-1))-exp(GAMMAtr))^2
                    - KAPPAi*(exp(I)*exp(dZtr)/exp(I(-1))-exp(GAMMAtr))
                      *exp(I)*exp(dZtr)/exp(I(-1)))
    + BETA*exp(Pk(+1))*exp(Zi(+1))
      *(1/exp(dZtr(+1))*(exp(C)-CHIc*exp(C(-1))/exp(dZtr))
        /(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1))))
      *KAPPAi*(exp(I(+1))*exp(dZtr(+1))/exp(I)-exp(GAMMAtr))
      *(exp(I(+1))*exp(dZtr(+1))/exp(I))^2;


//////////////////////////////// Production ///////////////////////////////
///////////////////////////////////////////////////////////////////////////
// MPK = MPL
exp(W)*(1+v*(exp(R)-1))/exp(Rk)
= (1-ALPHAk)/ALPHAk*exp(K(-1))*exp(ut)/exp(N)/exp(dZtr);

// MC
exp(MC) = 1/exp(Za)*(exp(W)*(1+v*(exp(R)-1)))^(1-ALPHAk)*exp(Rk)^ALPHAk
          *(1-ALPHAk)^(ALPHAk-1)*ALPHAk^(-ALPHAk);

exp(Y) = exp(Za)*(exp(K(-1))*exp(ut))^ALPHAk
         *exp(N)^(1-ALPHAk)*exp(dZtr)^(-ALPHAk)/exp(SS);

exp(Y) = exp(Ch) + exp(Ih) + exp(Xh) + exp(G)
         + KAPPAi/2*(exp(I)*exp(dZtr)/exp(I(-1))-exp(GAMMAtr))^2*exp(I)
         + 1/exp(Zi)*(GAMMA1*(exp(ut)-1)+GAMMA2*(exp(ut)-1)^2)*exp(ut)*exp(K(-1));

exp(G) = Gss*exp((Y))*exp(Zg)*(exp(Y)/exp(steady_state(Y)))^(-PHIg);


//////////////////////////////  UIP condition /////////////////////////////
///////////////////////////////////////////////////////////////////////////
exp(R) = exp(Rf)*exp(Zq)*exp(S(+1))/(exp(S));
exp(Rf) = exp(-UIPf*(exp(S)*exp(F)-(exp(steady_state(Y))*UIPy))
              -UIPr*(exp(Rstar)-exp(steady_state(Rstar))-(exp(R)-exp(steady_state(R)))))
          *exp(Zcp)*(exp(Rstar));   //exp(Rstar)-exp(steady_state(Rstar)) instead of Rf
//        *exp(Zcp)*exp(Rw)*(exp(Rstar));
// exp(Rstar)-exp(steady_state(Rstar)) instead of exp(Rf)-exp(steady_state(Rf))

exp(Rw) = (1+GAMMAtr)/BETA;
exp(Pxf)  = exp(PPx);//*exp(S));                                         
exp(Pm)  = exp(PPm)*(exp(PIstar))*exp(S);  //??
//exp(Pm)  = exp(Pm(-1))*(exp(PIstar));    //??
exp(Q) = exp(PPm)*(exp(PIstar))*exp(S)/exp(P);                          //??
//exp(Q) = exp(Q(-1))*(exp(PIstar))*exp(S)/exp(S(-1))*exp(P(-1))/exp(P);    //
//exp(Q) = exp(Rf)*exp(Q(+1))*exp(dP(+1))/exp(PIstar(+1))/exp(R);           //
//exp(Q) = exp(S)*exp(PIstar)/exp(P);
exp(PPx) = 1;
exp(PPm) = 1;
exp(M) = exp(Cf) + exp(If) + exp(Xf);


///////////////////////////  Taylor rule ////////////////////////
/////////////////////////////////////////////////////////////////
log(exp(R)) = er + TAYLORr*log(exp(R(-1)))
              + (1-TAYLORr)*(log(steady_state(exp(R)))+(1+TAYLORpi)*log(exp(dP))
                             + TAYLORy*log(exp(Y)/steady_state(exp(Y))));
                             

                             

/////////////////////////////////// shock processes //////////////////////
// permanent tech. shock: ln(A_t) = ln(A_t-1) + GAMMAtr*exp(Zg) + ea //
dZtr    = RHOtr *dZtr(-1)   + (1-RHOtr)*GAMMAtr + etr;                                                
Za      = RHOa  *Za(-1)     + ea;
Zi      = RHOi  *Zi(-1)     + ei;
Zc      = RHOc  *Zc(-1)     + ec;
Zpi     = RHOpi *Zpi(-1)    + epi;
Znw     = RHOnw *Znw(-1)    + enw;
Zpif    = RHOpif*Zpif(-1)   + epif;
Zx      = RHOx  *Zx(-1)     + ex;
Zg      = RHOg  *Zg(-1)     + eg;
Zw      = RHOw  *Zw(-1)     + ew;
Zcp     = RHOcp *Zcp(-1)    + ecp;
//Zmc     = (1-RHOmc)*PSIh    + RHOmc   *Zmc(-1)   + emc;
Zmc     = RHOmc *Zmc(-1)    + emc;
Zq      = RHOq *Zq(-1)    + eq;


////////////////////////////// Foreign VAR/////////////////////////////////
Ystar   = a11*Ystar(-1)  + a12*Rstar(-1) + a13*PIstar(-1) + eystar;
Rstar   = a21*Ystar(-1)  + a22*Rstar(-1) + a23*PIstar(-1) + erstar;
PIstar  = a31*Ystar(-1)  + a32*Rstar(-1) + a33*PIstar(-1) + epistar;

 

///* fb   ////// Flexible Price block /////////////////////////
//(1) c, r = R/dP(+1)
exp(dZtr(+1))*exp(Zc)/(exp(c)-CHIc*exp(c(-1))/exp(dZtr))
= BETA*exp(Zc(+1))/(exp(c(+1))-CHIc*exp(c)/exp(dZtr(+1)))*exp(r);

//(2) w = W/P
exp(w) = PSIw/(PSIw-1)*exp(Zw)*(exp(c)-CHIc*exp(c(-1))/exp(dZtr))*exp(n)^ETA;
//exp(w) = (exp(c)-CHIc*exp(c(-1))/exp(dZtr))*exp(n)^ETA;

//(3) ch, cf, ph, pcf
exp(ch) = ALPHAc*exp(ph)^(-XIc)*exp(c);

exp(cf) = (1-ALPHAc)*exp(pcf)^(-XIc)*exp(c);

1 = ALPHAc*(exp(ph))^(1-XIc) + (1-ALPHAc)*(exp(pcf))^(1-XIc);  // no shock on dP

//exp(ph) = PSIh/(PSIh-1)*exp(Zmc)*exp(mc);
exp(ph) = PSIh/(PSIh-1)*exp(mc);

exp(pcf) = PSIf/(PSIf-1)*exp(mcf); 

//(4) ih, ih, pi, pif
exp(ih) = ALPHAi*(exp(ph)/exp(pi))^(-XIi)*exp(i);

exp(if) = (1-ALPHAi)*(exp(pif)/exp(pi))^(-XIi)*exp(i);

exp(pi) = (ALPHAi*(exp(ph))^(1-XIi)+(1-ALPHAi)*(exp(pif))^(1-XIi))^(1/(1-XIi));

exp(pif) = PSIf/(PSIf-1)*exp(mcf);
//exp(pif) = exp(mcf);

//(5) x, xh, xh, pxf, px, f
exp(x) = exp(Zx)*(exp(x(-1))/exp(dZtr))^PHIx
         *(((1+v*(exp(r)-1))*exp(px)/exp(q))^(-XIh)*(exp(Ystar)))^(1-PHIx);
//exp(x) = (exp(x(-1))/exp(dZtr))^PHIx
//         *(((1+v*(exp(r)-1))*exp(px)/exp(q))^(-XIh)*(exp(Ystar)))^(1-PHIx);


exp(xh) = ALPHAx*(exp(ph)/exp(px))^(-XIx)*exp(x);

exp(xf) = (1-ALPHAx)*(exp(pxf)/exp(px))^(-XIx)*exp(x);

exp(px) = (ALPHAx*exp(ph)^(1-XIx)+(1-ALPHAx)*exp(pxf)^(1-XIx))^(1/(1-XIx));

exp(pxf) = exp(mcf);

exp(px)*exp(x) + exp(q)*exp(rf(-1))*exp(f(-1))/exp(dZtr)/exp(PIstar)
= exp(q)*exp(f) + exp(q)*(1+v*(exp(Rstar)-1))*(exp(cf) + exp(if) + exp(xf));

//(6) re, nw, k, Pk, rk, u, i
exp(re) = (exp(pk)*(1-DELTA) + exp(rk)*exp(u)
             - 1/exp(Zi)*(GAMMA1*(exp(u)-1)+GAMMA2/2*(exp(u)-1)^2))
          /exp(pk(-1));

exp(nw) = exp(Znw)*SP 
          *(exp(re)*exp(pk(-1))*exp(k(-1))/exp(dZtr)
            - (exp(nw(-1))/(exp(k(-1))*exp(pk(-1))))^(-KAPPArp)*exp(r(-1))
              *(exp(pk(-1))*exp(k(-1))/exp(dZtr)
                - exp(nw(-1))/exp(dZtr)));

exp(re(+1)) = (exp(nw)/(exp(k)*exp(pk)))^(-KAPPArp)*exp(r);

exp(k)-(1-DELTA)*exp(k(-1))/exp(dZtr)
= exp(Zi)*exp(i)*(1-KAPPAi/2*(exp(i)*exp(dZtr)/exp(i(-1))-exp(GAMMAtr))^2);

exp(rk) = (GAMMA1+GAMMA2*(exp(u)-1))/exp(Zi);

1 = exp(pk)*exp(Zi)*(1 - KAPPAi/2*(exp(i)*exp(dZtr)/exp(i(-1))-exp(GAMMAtr))^2
                      - KAPPAi*(exp(i)*exp(dZtr)/exp(i(-1))-exp(GAMMAtr))
                        *exp(i)*exp(dZtr)/exp(i(-1)))
    + BETA*exp(pk(+1))*exp(Zi(+1))
      *(1/exp(dZtr(+1))*(exp(c)-CHIc*exp(c(-1))/exp(dZtr))
        /(exp(c(+1))-CHIc*exp(c)/exp(dZtr(+1))))
      *KAPPAi*(exp(i(+1))*exp(dZtr(+1))/exp(i)-exp(GAMMAtr))
      *(exp(i(+1))*exp(dZtr(+1))/exp(i))^2;

//(7) w, rk, k, n, mc, mcf, y, ch, ih, xh, g
exp(w)*(1+v*(exp(r)-1))/exp(rk)
= (1-ALPHAk)/ALPHAk*exp(k(-1))*exp(u)/exp(n)/exp(dZtr);

exp(mc) = 1/exp(Za)*(exp(w)*(1+v*(exp(r)-1))/(1-ALPHAk))^(1-ALPHAk)
          *(exp(rk)/ALPHAk)^ALPHAk;

exp(mcf) = exp(q)*(1+v*(exp(Rstar)-1));

exp(y) = exp(Za)*(exp(k(-1))*exp(u))^ALPHAk*exp(n)^(1-ALPHAk)*exp(dZtr)^(-ALPHAk);

exp(y) = exp(ch) + exp(ih) + exp(xh) + exp(g)
         + KAPPAi/2*(exp(i)*exp(dZtr)/exp(i(-1))-exp(GAMMAtr))^2*exp(i)
         + 1/exp(Zi)*(GAMMA1*(exp(u)-1)+GAMMA2*(exp(u)-1)^2)
           *exp(u)*exp(k(-1));

exp(g) = Gss*exp((y))*exp(Zg)*(exp(y)/exp(steady_state(y)))^(-PHIg);

//(8) r, rf, f, Q, M, Rw
exp(r) = exp(rf)*exp(Zq)/exp(PIstar)*exp(q(+1))/(exp(q));
//exp(r) = exp(rf)*exp(q(+1))/(exp(q));
//exp(r) = exp(Rstar)/exp(PIstar);

exp(rf) = exp(- UIPf*(exp(q)*exp(f)-(exp(steady_state(y))*UIPy))
              - UIPr*(exp(Rstar)-exp(steady_state(Rstar))
                      - (exp(r)-exp(steady_state(r)))))
         *exp(Zcp)*(exp(Rstar));
          
        

exp(m) = exp(cf) + exp(if) + exp(xf);

//(9) GAP
Ygap = exp(Y) - exp(y);
Rgap = exp(R)/exp(dP(+1)) - exp(r);

//fe

// fcnb   ////// Flexible-Closed-No FA block: 
// variables ct, rt, wt, kt, nt, mct, yt, rkt, it, 
// financial accelerator variables: ret, nwt, pkt, utt

exp(dZtr(+1))/(exp(ct)-CHIc*exp(ct(-1))/exp(dZtr))
= BETA/(exp(ct(+1))-CHIc*exp(ct)/exp(dZtr(+1)))*exp(rt);

exp(wt) = PSIw/(PSIw-1)*(exp(ct)-CHIc*exp(ct(-1))/exp(dZtr))*exp(nt)^ETA;

1/exp(mct) = PSIh/(PSIh-1);

// 2 equations are modified and 4 equations are added for financial accelerator
// FA block starts
//exp(kt) - (1-DELTA)*exp(kt(-1))/exp(dZtr) = exp(it);
exp(kt)-(1-DELTA)*exp(kt(-1))/exp(dZtr)
= exp(Zi)*exp(it)*(1-KAPPAi/2*(exp(it)*exp(dZtr)/exp(it(-1))-exp(GAMMAtr))^2);

// exp(rt) = exp(rkt) + (1-DELTA);
exp(rkt) = (GAMMA1+GAMMA2*(exp(utt)-1))/exp(Zi);

exp(ret) = (exp(pkt)*(1-DELTA) + exp(rkt)*exp(utt)
             - 1/exp(Zi)*(GAMMA1*(exp(utt)-1)+GAMMA2/2*(exp(utt)-1)^2))
          /exp(pkt(-1));

exp(nwt) = exp(Znw)*SP 
          *(exp(ret)*exp(pkt(-1))*exp(kt(-1))/exp(dZtr)
            - (exp(nwt(-1))/(exp(kt(-1))*exp(pkt(-1))))^(-KAPPArp)*exp(rt(-1))
              *(exp(pkt(-1))*exp(kt(-1))/exp(dZtr)
                - exp(nwt(-1))/exp(dZtr)));

exp(ret(+1)) = (exp(nwt)/(exp(kt)*exp(pkt)))^(-KAPPArp)*exp(rt);

1 = exp(pkt)*exp(Zi)*(1 - KAPPAi/2*(exp(it)*exp(dZtr)/exp(it(-1))-exp(GAMMAtr))^2
                      - KAPPAi*(exp(it)*exp(dZtr)/exp(it(-1))-exp(GAMMAtr))
                        *exp(it)*exp(dZtr)/exp(it(-1)))
    + BETA*exp(pkt(+1))*exp(Zi(+1))
      *(1/exp(dZtr(+1))*(exp(ct)-CHIc*exp(ct(-1))/exp(dZtr))
        /(exp(ct(+1))-CHIc*exp(ct)/exp(dZtr(+1))))
      *KAPPAi*(exp(it(+1))*exp(dZtr(+1))/exp(it)-exp(GAMMAtr))
      *(exp(it(+1))*exp(dZtr(+1))/exp(it))^2;
// FA block ends


exp(wt) = exp(rkt)*(1-ALPHAk)/ALPHAk*exp(kt(-1))/exp(nt)/exp(dZtr);

exp(mct) = (exp(wt)/(1-ALPHAk))^(1-ALPHAk)*(exp(rkt)/ALPHAk)^ALPHAk;

exp(yt) = exp(kt(-1))^ALPHAk*exp(nt)^(1-ALPHAk)*exp(dZtr)^(-ALPHAk);

exp(yt) = exp(ct) + exp(it);

Ygapt = exp(Y) - exp(yt);

Rgapt = exp(R)/exp(dP(+1)) - exp(rt);

//(exp(rk) + (1-DELTA))*exp(PIstar)*exp(q_f) = exp(rf)*exp(q_f(+1));
//exp(rf) = (exp(Rstar));
//exp(q_f)*exp(rf(-1))*exp(f(-1)) = exp(q_f)*exp(f)*exp(dZtr)*exp(PIstar);


////fcne

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////


///////////////////////// Measurement Equations ////////////////
Y_obs       =   100*(Y-Y(-1) + dZtr) - 100*GAMMAtr;
R_obs       =   400*(R - steady_state(R));                      
C_obs       =   100*(C-C(-1) + dZtr) - 100*GAMMAtr;
dP_obs      =   (dP - steady_state(dP)); 
Inv_obs     =   100*(I-I(-1) + dZtr) - 100*GAMMAtr;
G_obs       =   100*(G-G(-1) + dZtr) - 100*GAMMAtr;
NN_obs      =   100*(NW-NW(-1));  
dPcf_obs    =   100*(dPcf) + 0.94;
//dPcf_obs    =   100*((dPcf)); 
dPx_obs     =   100*(dPx) - 0.258; 
//dPx_obs     =   100*((dPx)); 
X_obs       =   100*(X-X(-1) + dZtr) - 100*GAMMAtr;
S_obs       =   100*(Q-Q(-1));
dPi_obs     =   400*(dPi-steady_state(dPi));
SP_obs      =   400*(Re-R - (steady_state(Re) - steady_state(R)))                           ;  
Ystar_obs   =   100*(Ystar-Ystar(-1));
Rstar_obs   =   400*(Rstar); // steady state value or Rstar_obs  = 2.88, based on data
Pistar_obs  =   400*(PIstar); // steady state value or Rstar_obs  = 2.17. based on data
//W_obs       =   100*(W-W(-1) + dZtr) - 100*GAMMAtr;
/////////////////////////////////////////////////////////////////////
dY = (Y-Y(-4))*100;

end;


observation_trends;
Y_obs (GAMMAtr);
C_obs (GAMMAtr);
//Inv_obs (GAMMAtr);
//X_obs (GAMMAtr);
//G_obs (GAMMAtr);
//W_obs (GAMMAtr);

end;


initval;
Y           = 	 0.572716;
C           =    0.0854136;
N           =	-0.184133;
R          	=	 log((1+GAMMAtr)/BETA);//0.0125045;
W          	=	 0.161546;
Rk         	=	-2.82762;
K          	=	 2.11885;
dZtr        =	 GAMMAtr;
dP         	=	 0;
G          	=	-1.03672;
MC         	=	-0.182322;
Phstar      =	 0;
X1         	=	 1.74009;
X2         	=	 1.92241;
SS         	=	 0;
Q        	=	-0.286734;
F          	=	-0.227666;
Rf         	=	 0.0125045;
S          	=	-0.286735;
Pxf         =	 1.7212e-06;
Pm         	=	-0.286734;
X          	=	-0.631171;
PPx        	=	 4.00498e-07;
PPm        	=	 4.00498e-07;
Rw        	=	 0.0124597;
Pk          =	 2.23411e-11;
I        	=	-1.25595;
Re         	=	 0.0348178;
NW         	=	 1.89572;
Ch         	=	-0.24304;
Cf         	=	-1.16408;
P         	=	-0.0351037;
Pcfstar      =	-0.155084;
dPh        	=	  0;
dPcf        =	-0.0990233;
X1f        	=	-0.97575;
X2f        	=	-0.638344;
Ih       	=	-2.65812;
If       	=	-1.5351;
Pi       	=	-0.06606;
Pifstar     =	-0.0919085;
dPi      	=	 0;
dPif    	=	-0.0432459;
X1if    	=	-0.0223043;
X2if    	=	-0.0223043;
M        	=	-0.639335;
Zmc        	=	 0;
Zq        	=	 0;

/*
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
*/

Wstar       =	 0.161546;
X1w        	=	 1.19985;
ut         	=	 0.0632151;
Xh          =    log(0.2);
Xf          =    log(0.2);
Px          =    log(1);


///*//fb   ////// Flexible Price block /////////////////////////
c           =   0.0854136;
r           =   log((1+GAMMAtr)/BETA);
w           =   0.161546;
n           =   -0.184133;
ch          =   -0.24304;
cf          =   -1.16408;
pcf         =   -0.155084; ////
ph          =   0; ////
mc          =   -0.182322;
q           =   -0.286734;
ih          =   -2.65812;
if          =   -1.5351;
pi          =   -0.06606;
i           =   -1.25595;
pif         =   -0.0919085;///////////
re          =   0.0348178;
pk          =   2.23411e-11;
rk          =	-2.82762;
u           =   0.0632151;
nw          =   1.89572;
k           =   2.11885;
y           =   0.572716;
xh          =   log(0.2);
g           =   -1.03672;
x           =   -0.631171;
px          =   log(1);
xf          =   log(0.2);
rf          =   0.0125045;
f           =   -0.227666;
m           =   -0.639335;
Ygap        =   0;
Rgap        =   0;
//fe

ct          =   1.02618; 
rt          =   log((1+GAMMAtr)/BETA);
wt          =   0.44562;
kt          =   3.41132;
nt          =   0.327734;
mct         =   -0.168351;
yt          =   1.34218;
rkt         =   -3.33665;
it          =   0.0363352;
Ygapt       =   -0.727714;
Rgapt       =   0;

end;

//steady;
steady(solve_algo=4);
resid;
check;

shocks;
//  var etr; stderr SIGMAtr;
//  var ea;  stderr SIGMAa;
//  var er; stderr SIGMAr;
//  var ei; stderr SIGMAi;
//  var epi; stderr SIGMApi;
//  var es; stderr SIGMAs;
  var ecp; stderr SIGMAcp;
//  var ew; stderr SIGMAw;

end;
stoch_simul(order=1, irf=40) Y, C, I, dP, R, X, S, F, Ygap, Rgap, y, r, i, f, yt, rt;

//write_latex_dynamic_model;


/////////////////////////////////////////
//      Bayesian Estimation            //
//    (2001. 1/4~ 2012. 3/4)           //
/////////////////////////////////////////  

varobs Y_obs, R_obs, dP_obs, C_obs, S_obs,
//Inv_obs, NN_obs, X_obs,  
Ystar_obs,Rstar_obs, Pistar_obs ; //  W_obs
// G_obs, dPx_obs;
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
    PHIx,       beta_pdf,       0.5,     0.1;
    XIh,        gamma_pdf,      0.5,     0.1;
    ALPHAc,     beta_pdf,       0.6,     0.1;
    ALPHAi,     beta_pdf,       0.6,     0.1;     
    XIc,        gamma_pdf,      2,       0.5;
    XIi,        gamma_pdf,      2,       0.5;

    ALPHAx,     beta_pdf,       0.5,     0.1;
    XIx,        gamma_pdf,      2.5,     0.5;
    GAMMA2,     beta_pdf,       0.5,     0.1;     
    UIPr,       beta_pdf,       0.5,     0.1;     
    PHIg,       gamma_pdf,      0.5,     0.1;
    
    RHOa,       beta_pdf,       0.8,    0.1;
    RHOc,       beta_pdf,       0.8,    0.1;
    RHOtr,      beta_pdf,       0.8,    0.1;
    RHOpi,      beta_pdf,       0.8,    0.1;
    RHOi,       beta_pdf,       0.8,    0.1;
    RHOnw,      beta_pdf,       0.8,    0.1;
    RHOpif,     beta_pdf,       0.8,    0.1;
    RHOx,       beta_pdf,       0.8,    0.1;
    RHOpii,     beta_pdf,       0.8,    0.1;
    RHOg,       beta_pdf,       0.8,    0.1;
    RHOw,       beta_pdf,       0.8,    0.1;
    RHOcp,      beta_pdf,       0.8,    0.1;
    RHOmc,      beta_pdf,       0.8,    0.1;
    RHOq,       beta_pdf,       0.8,    0.1;


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
    stderr etr,     inv_gamma_pdf, 0.01,  2;
    stderr epi,     inv_gamma_pdf, 0.01,  2;
    stderr ei,      inv_gamma_pdf, 0.01,  2;
    stderr enw,     inv_gamma_pdf, 0.01,  2;
    stderr epif,    inv_gamma_pdf, 0.01,  2;
    stderr ex,      inv_gamma_pdf, 0.01,  2;
    stderr eg,      inv_gamma_pdf, 0.01,  2;
    stderr emc,     inv_gamma_pdf, 0.01,  2;
    stderr eystar,  inv_gamma_pdf, 0.01,  2;
    stderr erstar,  inv_gamma_pdf, 0.01,  2;
    stderr epistar, inv_gamma_pdf, 0.01,  2;
    stderr ew,      inv_gamma_pdf, 0.01,  2;
    stderr ecp,     inv_gamma_pdf, 0.01,  2;
    stderr epx,     inv_gamma_pdf, 0.01,  2;
    stderr eq,      inv_gamma_pdf, 0.01,  2;

end;


//be
estimation(
datafile     = data_201444_cpi, //data_201344_sa_soe_dmean.xls, //,
first_obs    = 1, 
mh_replic    = 2000,
//mode_file    = BOKDSGE2014_0623_TR_UIP_W_markup_mode,
mh_nblocks   = 1,
//forecast     = 12, 
//nobs       = [5:50],
mh_drop      = 0.4, 
mh_jscale    = 0.20, 
mode_compute = 4,
//mode_check, 
plot_priors  = 0,
smoother, 
//endogenous_prior,
//,mode_check
//filtered_vars,
nograph) Y_obs, R_obs,  dP_obs, C_obs, Inv_obs, NN_obs, dPcf_obs, X_obs, dPi_obs, 
S_obs, Ystar_obs, Rstar_obs, Pistar_obs dZtr G_obs, Za,  Zi, Zc,  
Zpi, Znw, Zpif, Zx, Zg, Zw, Zcp, Zmc, Ygap, Rgap, y, r, yt, rt, Ygapt, Rgapt;

shock_decomposition   
Y_obs, C_obs, R_obs, dP_obs, S_obs, Ygap, Rgap, y, r, yt, rt, Ygapt, 
Rgapt;
//be