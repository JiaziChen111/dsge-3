%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  New Keynesian model with Unit root shock process  %%
%%  2014.4.16                                        %%
%%  Written by Byoungho Bae                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

var   Y, C, N, R, W, Rk, K, dA, dP, G                                                                                                                                   // detrended variable (cf: Y = Y_trend/Z_t  where Z_t - Z_t-1 + gamma + ez)
      MC, P, X1, X2, SS, ERE                                                                                                                                            // variables required for price stickiness 
      F, RR, S, Px, Pm, X, PPx, PPm, RRw                                                                                                                                // variables for open economy
      Q, Inv, Re, NN,                                                                                                                                                   // variables for Financial Accelerator
      Cd, Cf, Pc, Pf, dPd, dPf, PPf, PP_invf, X1f, X2f, Invd, Invf, Pinv, P_invf, dPinv, dP_invf, X1_invf, X2_invf, Imp                                                 // variables for composite goods
      Ystar, Rstar, Pistar                                                                                                                                              // foreign VAR
      Za,  Zinv, Zc,  Zinf, Zn, Zinff, Zx, Zg, Zw, Zcp, Zmc,                                                                                                            // shock process 
      Y_obs, R_obs,  dP_obs, C_obs, Inv_obs, NN_obs, dPf_obs, X_obs, dPinv_obs, S_obs, Ystar_obs, Rstar_obs, Pistar_obs, W_obs, G_obs, SP_obs, dPx_obs                  // Observation variables
       GAPeff, Y_e, Ce, Ne, LAMe, We, R_e,Ke, Inve,  Rke, Qe                                                                                                            // variables for RBC block (for gap estimation)
      w, X1w, X1ww,X2w,X2ww, ut, Xd, Xm, Pxx, dPx 
      Ystar_or;                                                                                                                                                         // variables for sticky wage and utilization    

varexo er, 
       ea, 
       etrend, 
       einv, 
       ec, 
       einf, 
       en, 
       einff, 
       ex, 
       eg, 
       epr, 
       eystar, 
       erstar, 
       epistar, 
       ecp, 
       ew
       emc epx; 

////////////// define parameters //////////////
parameters bet, del, psi, alpk,  chic, eps, epsf, the, thef, rhop, rhoy, rhor, phi, zet, psi_x, psi_px, chik, gamma, gamma1, omega, 
           alp_c, rho_c, alp_i, rho_i, eta, del_e, iota, pi_target, thew, epsw, v, gam, gam2, zet2, alpx, etax
           rhogam, rhoa, rhoinv, rhoc,  rhoinf, rhon, rhoinff, rhox, rhopinv, g_ss, rhog, rhopr, rhow, rhocp, rhomc, 
           SIGMAr, SIGMAa, SIGMAinv, SIGMAc,  SIGMAinf, SIGMAn, SIGMAinff, SIGMAx, SIGMAg, SIGMAw,  SIGMAcp,  
           a11, a12, a13, a21, a22, a23, a31, a32, a33 phi_g;

bet     = 0.999;                // 할인률
del     = 0.025;                // 자본 감가상각률
gamma   = 0.0095;               // 영구적 기술충격내 성장률
chic    = 0.6909;               // 소비습관(consumption habit) 계수
chik    = 3.339;                // 투자조정 계수
del_e   = 0.025;       
psi     = 2;                    // parameter on labor function
alpk    = 0.33;                 // 생산함수내 자본투여 비율
the     = 0.807;                // 국내재 가격 조정 확률 (1- theta)
thef    = 0.546;                // 수입재국내재 가격 조정 확률 (1- theta)
thew    = 0.541;                // 임금 조정확률
rhop    = 0.4019;               // 테일러 준칙내 인플레이션 반응 계수
rhoy    = 0.312;                // 테일러 준칙내 생산갭 반응 계수
rhor    = 0.9407;               // 테일러 준칙내 금리 지속성 계수
eps     = 6.454;                // 국내재 대체탄력성
epsf    = 5.175;                // 수입재 대체탄력성
epsw    = 6.091;                // 가계 노동공급의 대체탄력성   
iota    = 0.3162;               // Price indexation 계수 
psi_x   = 0.4880;               // 수출함수 지속성 계수
psi_px  = 0.3945;               // 수출함수내 환율 탄력성 계수

phi     = 0.165;                 // UIP 조건내 계수
zet     = 0.03;                 // UIP 조건내 해외부채에 대한 탄력성 계수
zet2    = 0.5;
eta     = 1.2;                  // 노동함수내 Frisch 계수
gamma1  = 0.975;                // 기업가 생존확률
omega   = 0.05;                 // 리스크 프리미엄 반응계수
alp_c   = 0.56723;               // 소비복합재 국내재 비율
rho_c   = 3.955;                // 복합소비재 대체 탄력성
alp_i   = 0.35088;               // 투자복합재 국내재 비율 
rho_i   = 1.673;                // 투자복합재 대체 탄력성
g_ss    = 0.20; 
pi_target = 1.005;              // 물가목표수준 (1.005^4 = 1.020)   
v       = 1;                    // 운전자금(working capital) 대출 비율
gam     = 0.02;                 // capital utilization 계수
gam2    = 1.2;                  // capital utilization 계수
alpx    = 0.45550;
etax    = 1.488;
phi_g   = 1;

rhogam  = 0.8180;                  // 충격구조 지속성 계수
rhoa    = 0.9459;
rhoinv  = 0.9273;
rhoc    = 0.6902;
rhoinf  = 0.3224;
rhon    = 0.7981;
rhoinff = 0.9325;
rhox    = 0.8402;
rhopinv = 0.6705;
rhog    = 0.8052;
rhow    = 0.4684;
rhocp   = 0.90716;
rhomc   = 0.0;

SIGMAr   = 0.0014;             // 충격 표준편차
SIGMAa   = 0.0091;    
SIGMAinv = 0.1186;
SIGMAc   = 0.0261;
SIGMAinf = 0.0122;
SIGMAn   = 0.0034;
SIGMAinff= 0.0502;
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
///////////////////////////// 소비결정식 /////////////////////////////////////
exp(dA(+1))* exp(Zc)*1/(exp(C)-chic*exp(C(-1))/exp(dA)) = bet*exp(Zc(+1))*1/(exp(C(+1))-chic*exp(C)/exp(dA(+1))) *exp(R)/exp(dP(+1));   


/////////////  임금경직성 관련 균형식 /////////////////  
exp(X1w) = exp(w)^(-epsw*(eta+1))*exp(X1ww);
exp(X1ww) = exp(W)^(epsw*(eta+1))*exp(N)^(eta+1) + bet*thew*exp(dA(+1))^(epsw*(eta+1))*exp(X1ww(+1));
exp(X2w) = exp(w)^(-epsw)*exp(X2ww);
exp(X2ww) = 1/(exp(C)-chic*exp(C(-1))/exp(dA))*exp(W)^epsw*exp(N) + bet*thew*exp(dA(+1))^(epsw-1)*exp(X2ww(+1));
exp(w)= epsw/(epsw-1)*exp(Zw)*exp(X1w) / exp(X2w);
exp(W)   = (thew*exp(W(-1)/exp(dA))^(1-epsw) + (1-thew)*exp(w)^(1-epsw))^(1/(1-epsw));                                                     

/////////////////////////// 소비복합재 관련 조건식 ///////////////////////////////////////////
exp(Cd) = alp_c*exp(Pc)^rho_c*exp(C);                        
exp(Cf) = ((1-alp_c)*(exp(Pc)/exp(PPf))^rho_c*exp(C));                
exp(Pc) = (alp_c + (1-alp_c)*(exp(PPf))^(1-rho_c))^(1/(1-rho_c));   
exp(dP)*exp(Zinf) = (exp(Pc)/exp(Pc(-1))*exp(dPd));                

///////////////////////////// 가격 경직성 조건식 (국내 및 수입 소비재) ////////////////////////////////////////////////////
exp(X1)         = exp(Y)*exp(MC) + the*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *exp(dPd(+1))^(Zmc+1)*exp(dPd)^((iota)*(-Zmc)) *pi_target^((-Zmc)*(1-iota))*exp(X1(+1))  ;
exp(X2)         = exp(Y)         + the*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *exp(dPd(+1))^(Zmc)  *exp(dPd)^((iota)*(1-Zmc))*pi_target^((1-Zmc)*(1-iota))*exp(X2(+1))  ;
exp(P)          = Zmc/(Zmc-1)* exp(X1)/ exp(X2);
1               = exp(dPd(-1))^(iota*(1-Zmc))*pi_target^((1-iota)*(1-Zmc))*the*exp(dPd)^(Zmc-1) + (1-the)*exp(P)^(1-Zmc);
exp(SS)         = (exp(dPd(-1))^iota*pi_target^(1-iota))^(-eps)*the*exp(dPd)^(eps)*exp(SS(-1)) + (1-the)*exp(P)^(-eps);

exp(X1f)         = exp(Cf)*exp(ERE)*exp(Pc)*(1+v*(exp(R)-1))*exp(Pistar) + thef*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *(exp(dPf(+1)))^(epsf)*(exp(dPf))^((iota)*(-epsf)) *pi_target^(-epsf*(1-iota))*exp(dPd(+1)) *exp(X1f(+1))  ;
exp(X2f)         = exp(Cf)                                               + thef*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *(exp(dPf(+1)))^(epsf)*(exp(dPf))^((iota)*(1-epsf)) *pi_target^((1-epsf)*(1-iota)) *exp(X2f(+1))  ;
exp(Pf)          = epsf/(epsf-1)* exp(X1f)/ exp(X2f);
1                = exp(dPf(-1))^(iota*(1-epsf))*pi_target^((1-iota)*(1-epsf))* thef*exp(dPf)^(epsf-1) + (1-thef)*(exp(Pf)/exp(PPf))^(1-epsf);
exp(dPf)*exp(Zinff)         = exp(PPf)/exp(PPf(-1))*exp(dPd);
 
/////////////////////////// 투자복합재 관련 조건식 ////////////////////////////////////////////
exp(Invd) = alp_i*exp(Pinv)^rho_i*exp(Inv);                        
exp(Invf) = (1-alp_i)*(exp(Pinv)/exp(PP_invf))^rho_i*exp(Inv);                
exp(Pinv) = (alp_i + (1-alp_i)*(exp(PP_invf))^(1-rho_i))^(1/(1-rho_i));   
exp(dPinv) = exp(Pinv)/exp(Pinv(-1))*exp(dPd);                

exp(X1_invf)    = exp(Invf)*exp(ERE)*exp(Pc)*(1+v*(exp(R)-1))*exp(Pistar) + thef*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *(exp(dP_invf(+1)))^(epsf) *(exp(dP_invf))^((iota)*(-epsf)) *pi_target^(-epsf*(1-iota))*exp(dPd(+1))*exp(X1_invf(+1))  ;
exp(X2_invf)    = exp(Invf)                                               + thef*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *(exp(dP_invf(+1)))^(epsf) *(exp(dP_invf))^((iota)*(1-epsf)) *pi_target^((1-epsf)*(1-iota))*exp(X2_invf(+1))  ;
exp(P_invf)     = epsf/(epsf-1)* exp(X1_invf)/ exp(X2_invf);
1               =  exp(dP_invf(-1))^(iota*(1-epsf))*pi_target^((1-iota)*(1-epsf))*thef*(exp(dP_invf))^(epsf-1) + (1-thef)*(exp(P_invf)/exp(PP_invf))^(1-epsf);
exp(dP_invf)*exp(Zinff)    = exp(PP_invf)/exp(PP_invf(-1))*exp(dPd);

//////////////////////////////////// 금융가속기 (Financial Accelerator) ////////////////////////////////////////////
exp(Re)          = (exp(Q)*(1-del) + exp(Rk)*exp(ut) - 1/exp(Zinv)*(gam*(exp(ut)-1) + gam2/2*(exp(ut)-1)^2) )*exp(dPd)/exp(Q(-1));
exp(NN)          = exp(Zn)*gamma1*( exp(Re)*exp(Q(-1))*exp(K(-1))/exp(dA) - ( exp(NN(-1))/ (exp(K(-1))* exp(Q(-1))) )^(-omega)*exp(R(-1))/exp(dP)*(exp(Q(-1))*exp(K(-1))/exp(dA)- exp(NN(-1))/exp(dA)));
exp(Re(+1))      = (exp(NN)/(exp(K)*exp(Q)))^(-omega)*exp(R)/exp(dP(+1)) ;

exp(Zinv)*exp(Inv) *(1-chik/2*(exp(Inv)*exp(dA)/exp(Inv(-1)) - exp(gamma))^2)  =  exp(K)- (1- del)*exp(K(-1))/exp(dA);                                     %% 투자 동태식 
exp(Rk)          = (gam + gam2*(exp(ut) -1))/ exp(Zinv);                                                                                                    %% 자본 가동률에 대한 일계조건식 Rk = A(ut)'
1                = exp(Q)*exp(Zinv)*(1-chik/2*(exp(Inv)*exp(dA)/exp(Inv(-1)) -exp(gamma))^2 -chik*(exp(Inv)*exp(dA)/exp(Inv(-1)) -exp(gamma))*exp(Inv)*exp(dA)/exp(Inv(-1)))
                 + bet*exp(Q(+1))*exp(Zinv(+1))*(1/exp(dA(+1))*(exp(C)-chic*exp(C(-1))/exp(dA))/(exp(C(+1))-chic*exp(C)/exp(dA(+1))))*chik*(exp(Inv(+1))*exp(dA(+1))/exp(Inv)-exp(gamma))*(exp(Inv(+1))*exp(dA(+1))/exp(Inv))^2;    %% 투자에 대한 일계조건식 (자본재 가격 결정식)


////////////////////////////////// 생산부문 //////////////////////////////////////////////////
exp(W)*(1+v*(exp(R)-1))/exp(Rk)        = (1-alpk)/alpk*exp(K(-1))*exp(ut)/exp(N)*1/exp(dA);                                                       // 자본 노동의 생산투입요소 비율
exp(MC)          = 1/exp(Za)*(exp(W)*(1+v*(exp(R)-1)))^(1-alpk)*exp(Rk)^alpk*(1-alpk)^(alpk-1)*alpk^(-alpk);                  // 마크업 정의식 

//exp(W)*(1+v*(exp(R)-1)) = exp(Za)*exp(MC)* (1-alpk)*(exp(K(-1))*exp(ut))^(alpk)   *exp(N)^(-alpk) *exp(dA)^(-alpk)/exp(SS);
//exp(Rk)                 = exp(Za)*exp(MC)* alpk*(exp(K(-1))*exp(ut))^(alpk-1)  *exp(N)^(1-alpk)*exp(dA)^(1-alpk)/exp(SS);
exp(Y)                  = exp(Za)*(exp(K(-1))*exp(ut))^alpk*exp(N)^(1-alpk)*exp(dA)^(-alpk)/exp(SS);

exp(Y)                  = exp(Cd) + exp(Invd) + exp(Xd) + chik/2*(exp(Inv)*exp(dA)/exp(Inv(-1))-exp(gamma))^2*exp(Inv) + exp(G) + 1/exp(Zinv)*(gam*(exp(ut)-1) + gam2*(exp(ut)-1)^2)*exp(ut)*exp(K(-1));
exp(G)                  = g_ss*exp((Y))*exp(Zg)*(exp(Y)/exp(steady_state(Y)))^(-phi_g);


/////////////////////////////////// 경제충격 프로세스 /////////////////////////////////////////////////////
dA         = rhogam*dA(-1) + (1-rhogam)*gamma + etrend;                                                // 영구적 기술충격 ln(A_t) = ln(A_t-1) + gamma*exp(Zg) + ea //
Za         = rhoa    *Za(-1)      + ea;
Zinv       = rhoinv  *Zinv(-1)    + einv;
Zc         = rhoc    *Zc(-1)      + ec;
Zinf       = rhoinf  *Zinf(-1)    + einf;
Zn         = rhon    *Zn(-1)      + en;
Zinff      = rhoinff *Zinff(-1)   + einff;
Zx         = rhox    *Zx(-1)      + ex;
Zg         = rhog    *Zg(-1)      + eg;
Zw         = rhow    *Zw(-1)      + ew;
Zcp        = rhocp   *Zcp(-1)     + ecp;
Zmc        = (1-rhomc)*eps + rhomc   *Zmc(-1)   + emc;


////////////////////////////////////// 통화준칙 (Taylor rule) /////////////////////////////////////////////////////////
log(exp(R))     =   rhor*log(exp(R(-1)))+ (1-rhor)*(log(steady_state(exp(R)))+ (1+rhop)*log(exp(dP))+ rhoy*log(exp(Y)/steady_state(exp(Y)))) +er ;


////////////////////////////////////// Foreign VAR/////////////////////////////////////////////////////////
Ystar           = a11*Ystar(-1)  + a12*Rstar(-1) + a13*Pistar(-1) + eystar;
Rstar           = a21*Ystar(-1)  + a22*Rstar(-1) + a23*Pistar(-1) + erstar;
Pistar          = a31*Ystar(-1)  + a32*Rstar(-1) + a33*Pistar(-1) + epistar;

 
//////////////////////////////////////// Export function & BOP  //////////////////////////////////////////
exp(X)          = exp(Zx)*(exp(X(-1))/exp(dA))^psi_x *(( (1+v*(exp(R)-1))*exp(Pxx)/exp(S))^(-psi_px)*(exp(Ystar)))^(1-psi_x);          // 수출함수
exp(Xd)         = alpx*exp(Pxx)^etax*exp(X);
exp(Xm)         = (1-alpx)*(exp(Pxx)/(exp(ERE)*exp(Pc)))^etax*exp(X);
exp(Pxx)        = (alpx + (1-alpx)*(exp(ERE)*exp(Pc))^(1-etax))^(1/(1-etax));
exp(dPx)*exp(epx) = exp(Pxx)/exp(Pxx(-1))*exp(dPd);
exp(ERE)*exp(Px)*exp(Pc)*exp(X) + exp(S)*exp(RR(-1))*exp(F(-1))*exp(Pc)/exp(dA)/exp(dP)  =  exp(S)*exp(F)*exp(Pc) + exp(ERE)*exp(Pc)*exp(Cf)*(1+v*(exp(R)-1)) + exp(ERE)*exp(Pc)*exp(Invf)*(1+v*(exp(R)-1)) + exp(ERE)*exp(Pc)*exp(Xm)*(1+v*(exp(R)-1));

////////////////////////////////////////  UIP condition //////////////////////////////////////////
exp(R)          = exp(RR)*exp(S(+1))/(exp(S));                                               // UIP condition
exp(RR)         = exp(RRw)*(exp(Rstar)) * exp(-zet*(exp(S)*exp(F)- (exp(steady_state(Y))*phi)) - zet2*(exp(RR)-exp(steady_state(RR))- (exp(R)-exp(steady_state(R)))))*exp(Zcp);                     // 

exp(RRw)        = (1+gamma)/bet;
exp(Px)         = exp(PPx);//*exp(S));                                         
exp(Pm)         = exp(PPm)*(exp(Pistar))*exp(S);                                            // 수입재 가격
exp(ERE)        = exp(PPm)*(exp(Pistar))*exp(S)/exp(Pc);                                    // 실질환율
exp(PPx)        = 1;
exp(PPm)        = 1;
exp(Imp)        = exp(Cf) + exp(Invf) + exp(Xm);
  


/////////////////////////////////////// Efficient Gap 추정을 위한 RBC block ////////////////////////////////////////////////////////
 exp(LAMe)                          = 1/(exp(Ce));                                                                                      // RBC (1)
 psi/(1-exp(Ne))                    = 1/exp(Ce)*exp(We); //psi*exp(Ne)^eta                    = exp(LAMe)*exp(We);                      // RBC (2)
 exp(dA(+1))* exp(Qe)*1/exp(Ce)     = bet* 1/exp(Ce(+1)) * (exp(Rke(+1)) +exp(Qe(+1))*(1-del_e) );                                      // RBC (3)
 exp(dA(+1))* 1/exp(Ce)             = bet*( 1/exp(Ce(+1)) *exp(R_e));                                                                   // RBC (4)
 exp(We)                            = (1-alpk)*exp(Za)*exp(Ke(-1))^(alpk)*exp(Ne)^(-alpk) *exp(dA)^(1-alpk);                            // RBC (5)
 exp(Rke)                           = alpk*exp(Za)*exp(Ke(-1))^(alpk-1)*exp(Ne)^(1-alpk)*exp(dA)^(1-alpk);                              // RBC (6)
 exp(Y_e)                           = exp(Za)*exp(Ke(-1))^alpk*exp(Ne)^(1-alpk)*exp(dA)^(-alpk);                                        // RBC (7)
 exp(Y_e)                           = exp(Ce) + exp(Ke) - (1-del)*exp(Ke(-1))*1/exp(dA)+ chik/2*(exp(Inve)/exp(Ke)-del)^2*exp(Ke);      // RBC (8)
 exp(Qe)                            = (1+chik*(exp(Inve)/ exp(Ke) - del_e));                                                            // RBC (9)
 exp(Inve)*exp(Zinv)                =  exp(Ke)- (1- del)*exp(Ke(-1))*1/exp(dA);                                                         // RBC (10)

 GAPeff = (exp(Y)-exp(steady_state(Y)))/exp(steady_state(Y))  - (exp(Y_e)-exp(steady_state(Y_e)))/exp(steady_state(Y_e));
 
 //////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////// 관측방정식(Measurement Equations) ////////////////
Y_obs       =   100*(Y-Y(-1) + dA)          - 100*gamma;
R_obs       =   400*(R - steady_state(R));                      
C_obs       =   100*(C-C(-1) + dA)          - 100*gamma;
dP_obs      =   100*(dP - steady_state(dP)); 
Inv_obs     =   100*(Inv-Inv(-1) + dA)      - 100*gamma;
G_obs       =   100*(G-G(-1) + dA)          - 100*gamma;
NN_obs      =   100*(NN-NN(-1));  
dPf_obs     =   100*(dPf) + 0.94;                 
dPx_obs     =   100*(dPx) - 0.258;                         
X_obs       =   100*(X-X(-1) + dA)          - 100*gamma;
S_obs       =   100*(ERE-ERE(-1)); 
dPinv_obs   =   400*(dPinv-steady_state(dPinv));
SP_obs      =   400*(Re-R -(steady_state(Re) - steady_state(R)))                           ;  
Ystar_obs   =   100*(Ystar-Ystar(-1));
Rstar_obs   =   400*(Rstar)                       ; // steady state value or Rstar_obs  = 2.88, based on data
Pistar_obs  =   400*(Pistar)                      ; // steady state value or Rstar_obs  = 2.17. based on data
W_obs       =   100*(W-W(-1)+ dA)           - 100*gamma;
/////////////////////////////////////////////////////////////////////

Ystar_or = (Y-Y(-4))*100;


end;

observation_trends;
Y_obs (gamma);
C_obs (gamma);
Inv_obs (gamma);
X_obs (gamma);
G_obs (gamma);
W_obs (gamma);

end;


initval;
Y           = 	 0.572716;
C           =    0.0854136;
N           =	-0.184133;
R          	=	 log((1+gamma)/bet);//0.0125045;
W          	=	 0.161546;
Rk         	=	-2.82762;
K          	=	 2.11885;
dA         	=	 gamma;
dP         	=	 0;
G          	=	-1.03672;
MC         	=	-0.182322;
P          	=	 0;
X1         	=	 1.74009;
X2         	=	 1.92241;
SS         	=	 0;
ERE        	=	-0.286734;
F          	=	-0.227666;
RR         	=	 0.0125045;
S          	=	-0.286735;
Px         	=	 1.7212e-06;
Pm         	=	-0.286734;
X          	=	-0.631171;
PPx        	=	 4.00498e-07;
PPm        	=	 4.00498e-07;
RRw        	=	 0.0124597;
Q          	=	 2.23411e-11;
Inv        	=	-1.25595;
Re         	=	 0.0348178;
NN         	=	 1.89572;
Cd         	=	-0.24304;
Cf         	=	-1.16408;
Pc         	=	-0.0351037;
Pf         	=	-0.155084;
dPd        	=	  0;
dPf        	=	-0.0990233;
X1f        	=	-0.97575;
X2f        	=	-0.638344;
Invd       	=	-2.65812;
Invf       	=	-1.5351;
Pinv       	=	-0.06606;
P_invf     	=	-0.0919085;
dPinv      	=	 0;
dP_invf    	=	-0.0432459;
X1_invf    	=	-0.0223043;
X2_invf    	=	-0.0223043;
Imp        	=	-0.639335;
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
w          	=	 0.161546;
X1w        	=	 1.19985;
ut         	=	 0.0632151;
Xd          =    log(0.2);
Xm          =    log(0.2);
Pxx         =    log(1);
end;

steady;
resid;
check;

shocks;
//  var etrend; stderr SIGMAtr;
//  var ea; stderr SIGMAa;
 // var er; stderr SIGMAr;
//  var einv; stderr SIGMAinv;
//  var einf; stderr SIGMAinf;
//  var es; stderr SIGMAs;
  var ecp; stderr SIGMAcp;
//  var ew; stderr SIGMAw;

end;
stoch_simul(order=1, irf=30) Y, C, Inv, dP, R, X,S, F;

//write_latex_dynamic_model;





/////////////////////////////////////////
//      Bayesian Estimation            //
//    (2001. 1/4~ 2012. 3/4)           //
/////////////////////////////////////////  

varobs Y_obs, R_obs, dP_obs, C_obs, Inv_obs, NN_obs, X_obs,S_obs, W_obs, Ystar_obs, Rstar_obs, Pistar_obs, G_obs, dPx_obs;//, dPf_obs;//,  G_obs;//;// dPf_obs;//, W_obs ;//dPf_obs,  dPinv_obs ;//,  Ystar_obs, Rstar_obs, Pistar_obs G_obs, W_obs, dPd_obs;//, N_obs;//, SP_obs;

estimated_params;
    eta,        gamma_pdf,      1.4,     0.2;
    eps,        gamma_pdf,      6,       1.0;
    epsf,       gamma_pdf,      6,       1.0;
    epsw,       gamma_pdf,      6,       1.0;
    the,        beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05 
    thef,       beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05
    thew,       beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05
    chic,       normal_pdf,     0.5,     0.1;
    chik,       normal_pdf,     4,       1.0;
    omega,      beta_pdf,       0.3,     0.05;    
    
    rhor,       beta_pdf,       0.8,     0.05;
    rhoy,       beta_pdf,       0.2,     0.05;
    rhop,       beta_pdf,       0.5,     0.1;   
    iota,       beta_pdf,       0.5,     0.1;
    psi_x,      beta_pdf,       0.5,     0.1;
    psi_px,     gamma_pdf,      0.5,     0.1;
    alp_c,      beta_pdf,       0.6,     0.1;
    alp_i,      beta_pdf,       0.6,     0.1;     
    rho_c,      gamma_pdf,      2,       0.5;
    rho_i,      gamma_pdf,      2,       0.5;

    alpx,       beta_pdf,       0.5,     0.1;
    etax,       gamma_pdf,      2.5,     0.5;
    gam2,       beta_pdf,       0.5,     0.1;     
    zet2,       beta_pdf,       0.5,     0.1;     
    phi_g,      gamma_pdf,      0.5,     0.1;
    
    rhoa,       beta_pdf,       0.8,    0.1;
    rhoc,       beta_pdf,       0.8,    0.1;
    rhogam,     beta_pdf,       0.8,    0.1;
    rhoinf,     beta_pdf,       0.8,    0.1;
    rhoinv,     beta_pdf,       0.8,    0.1;
    rhon,       beta_pdf,       0.8,    0.1;
    rhoinff,    beta_pdf,       0.8,    0.1;
    rhox,       beta_pdf,       0.8,    0.1;
    rhopinv,    beta_pdf,       0.8,    0.1;
    rhog,       beta_pdf,       0.8,    0.1;
    rhow,       beta_pdf,       0.8,    0.1;
    rhocp,      beta_pdf,       0.8,    0.1;

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
    stderr etrend,  inv_gamma_pdf, 0.01,  2;
    stderr einf,    inv_gamma_pdf, 0.01,  2;
    stderr einv,    inv_gamma_pdf, 0.01,  2;
    stderr en,      inv_gamma_pdf, 0.01,  2;
    stderr einff,   inv_gamma_pdf, 0.01,  2;
    stderr ex,      inv_gamma_pdf, 0.01,  2;
    stderr eg,      inv_gamma_pdf, 0.01,  2;
    stderr emc,     inv_gamma_pdf, 0.01,  2;
    stderr eystar,  inv_gamma_pdf, 0.01,  2;
    stderr erstar,  inv_gamma_pdf, 0.01,  2;
    stderr epistar, inv_gamma_pdf, 0.01,  2;
    stderr ew,      inv_gamma_pdf, 0.01,  2;
    stderr ecp,     inv_gamma_pdf, 0.01,  2;
    stderr epx,    inv_gamma_pdf, 0.01,  2;



end;

/*
estimation(
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
nograph ) Y_obs, R_obs,  dP_obs, C_obs, Inv_obs, NN_obs, dPf_obs, X_obs, dPinv_obs, S_obs, Ystar_obs, Rstar_obs, Pistar_obs dA G_obs W_obs Za,  Zinv, Zc,  Zinf, Zn, Zinff, Zx, Zg, Zw, Zcp;

shock_decomposition   Y_obs, C_obs, R_obs, dP_obs, Inv_obs, GAPeff;
     
*/
            
     