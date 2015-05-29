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
           a11, a12, a13, a21, a22, a23, a31, a32, a33, phi_g,



// steady_state variables //         
           Y_ss, C_ss, N_ss, R_ss, W_ss, Rk_ss, K_ss, dA_ss, dP_ss, G_ss, MC_ss, P_ss, X1_ss, X2_ss,
           SS_ss, ERE_ss, F_ss, RR_ss, S_ss, Px_ss, Pm_ss, X_ss, PPx_ss, PPm_ss, RRw_ss, Q_ss, Inv_ss,
           Re_ss, NN_ss, Cd_ss, Cf_ss, Pc_ss, Pf_ss, dPd_ss, dPf_ss, PPf_ss, PP_invf_ss, X1f_ss, X2f_ss,
           Invd_ss, Invf_ss, Pinv_ss, P_invf_ss, dPinv_ss, dP_invf_ss, X1_invf_ss, X2_invf_ss, Imp_ss, 
           Ystar_ss, Rstar_ss, Pistar_ss, Za_ss, Zinv_ss, Zc_ss, Zinf_ss, Zx_ss, Zg_ss, Zw_ss, Zcp_ss, Zn_ss
           Zmc_ss, Y_obs_ss, R_obs_ss, dP_obs_ss, C_obs_ss, Inv_obs_ss, NN_obs_ss, dPf_obs_ss, X_obs_ss,
           dPinv_obs_ss, S_obs_ss, Ystar_obs_ss, Rstar_obs_ss, Pistar_obs_ss, W_obs_ss, G_obs_ss, SP_obs_ss,
           dPx_obs_ss, GAPeff_ss, Y_e_ss, Ce_ss, Ne_ss, LAMe_ss, We_ss, R_e_ss, Ke_ss, Inve_ss, Rke_ss, Qe_ss,
           w_ss, X1w_ss, X1ww_ss, X2w_ss, X2ww_ss, ut_ss, Xd_ss, Xm_ss, Pxx_ss, dPx_ss, Ystar_or_ss 



;

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



//steady_state//

Y_ss = 1.13119084748 ;
C_ss = 0.60625951214 ;
N_ss = 0.37369265825 ;
R_ss = 0.01050050033 ;
W_ss = 0.17811217949 ;
Rk_ss = -2.81721866732 ;
K_ss =  2.64823304686 ;
dA_ss = 0.00950000000 ;
dP_ss = 0.00000000000 ;
G_ss = -0.47824706496 ;
MC_ss = -0.16913302902 ;
P_ss = -0.01361141238 ;
X1_ss = 2.51625764962 ;
X2_ss = 2.69821987121 ;
SS_ss = 0.00072508622 ;
ERE_ss = 0.04766369226 ;
F_ss = -0.84031646611 ;
RR_ss = 0.01050050033 ;
S_ss = 0.16690158463 ;
Px_ss = 0.00000172120 ;
Pm_ss = 0.16690198542 ;
X_ss = 0.01746058156 ;
PPx_ss = 0.00000040050 ;
PPm_ss = 0.00000040050 ;
RRw_ss = 0.01045766585 ;
Q_ss = 0.00000000002 ;
Inv_ss = -0.72675164538 ;
Re_ss = 0.03481780798 ;
NN_ss = 2.16188689389 ;
Cd_ss = 0.51085654666 ;
Cf_ss = -1.31090409273 ;
Pc_ss = 0.11923829287 ;
Pf_ss = 0.38817464729 ;
dPd_ss = 0.00000000000 ;
dPf_ss = 0.00000000008 ;
PPf_ss = 0.39221294140 ;
PP_invf_ss = 0.39221240374 ;
X1f_ss = -0.36582136226 ;
X2f_ss = -0.53927102871 ;
Invd_ss = -1.36832585679 ;
Invf_ss = -1.40932389346 ;
Pinv_ss = 0.24252049187 ;
P_invf_ss = 0.38817446823 ;
dPinv_ss = 0.00000000000 ;
dP_invf_ss = 0.00000044990 ;
X1_invf_ss = -0.46423932484 ;
X2_invf_ss = -0.63768881222 ;
Imp_ss = 0.00744009976 ;
Ystar_ss = 0.00000000000 ;
Rstar_ss = 0.00000000000 ;
Pistar_ss = 0.00000000000 ;
Za_ss = 0.00000000000 ;
Zinv_ss = 0.00000000000 ;
Zc_ss = 0.00000000000 ;
Zinf_ss = 0.00000000000 ;
Zn_ss = 0.00000000000 ; 
Zinff_ss = 0.00000000000 ;
Zx_ss = 0.00000000000 ;
Zg_ss = 0.00000000000 ;
Zw_ss = 0.00000000000 ;
Zcp_ss = 0.00000000000 ;
Zmc_ss = 6.45400000000 ;
Y_obs_ss = 0.00000000000 ;
R_obs_ss = 0.00000000000 ;
dP_obs_ss = 0.00000000000 ;
C_obs_ss = 0.00000000000 ;
Inv_obs_ss = 0.00000000000 ;
NN_obs_ss = 0.00000000000 ;
dPf_obs_ss = 0.94000000822 ;
X_obs_ss = 0.00000000000 ;
dPinv_obs_ss = 0.00000000000 ;
S_obs_ss = 0.00000000000 ;
Ystar_obs_ss = 0.00000000000 ;
Rstar_obs_ss = 0.00000000000 ;
Pistar_obs_ss = 0.00000000000 ;
W_obs_ss = 0.00000000000 ;
G_obs_ss = 0.00000000000 ;
SP_obs_ss = 0.00000000000 ;
dPx_obs_ss = -0.25800000000 ;
GAPeff_ss = 0.00000000000 ;
Y_e_ss = -0.02722764861 ;
Ce_ss = -0.40167832198 ;
Ne_ss = -1.10965997753 ;
LAMe_ss = 0.40167834884 ;
We_ss = 0.69145476234 ;
R_e_ss = 0.01050251301 ;
Ke_ss = 2.17993798893 ;
Inve_ss = -1.19504670042 ;
Rke_ss = -3.30632826186 ;
Qe_ss = 0.03031681247 ;
w_ss = 0.18011580569 ;
X1w_ss = 1.74675667913 ;
X1ww_ss = 4.16034449855 ;
X2w_ss = 1.74597887036 ;
X2ww_ss = 2.84306424283 ;
ut_ss = 0.03260590061 ;
Xd_ss = -0.63618613273 ;
Xm_ss = -0.70606406092 ;
Pxx_ss = 0.08918874193 ;
dPx_ss = 0.00000000000 ;
Ystar_or_ss = 0.00000000000 ;





/////////////////////////////////////////////
//               Model (F.O.C)             //
/////////////////////////////////////////////

model;
///////////////////////////// 소비결정식 /////////////////////////////////////

// %% eq (1)  
(Zc - Zc_ss) + (dA(+1) - dA_ss) + (-(C - C_ss) + (chic/exp(dA_ss))*(C(-1) - C_ss) - (chic/exp(dA_ss))*(dA - dA_ss))/(1 - (chic/exp(dA_ss))) = (-(C(+1) - C_ss) + (chic/exp(dA_ss))*(C - C_ss) - (chic/exp(dA_ss))*(dA(+1) - dA_ss))/(1 - (chic/exp(dA_ss))) + (R - R_ss)- (dP(+1) - dP_ss) + (Zc(+1) - Zc_ss)  ;


/////////////  임금경직성 관련 균형식 /////////////////  


// %% eq(2)
X1w = X1ww - epsw*(eta+1)*w; 

// %% eq(3)                                                                                                                                                                                              
((exp(W_ss)^(epsw*(eta+1)))*(exp(N_ss)^(eta+1)) + bet*thew*(exp(dA_ss)^(epsw*(eta+1)))*exp(X1ww_ss))*(X1ww - X1ww_ss) = (eta+1)*(exp(W_ss)^(epsw*(eta+1)))*(exp(N_ss)^(eta+1))*(epsw*(W - W_ss) + (N - N_ss)) + bet*thew*(exp(dA_ss)^(epsw*(eta+1)))*exp(X1ww_ss)*((epsw*(eta+1))*(dA(+1) - dA_ss) + (X1ww(+1) - X1ww_ss)) ;

// %% eq(4)                                                                                                                                                                                                                                                                                                                            
X2w = X2ww - epsw*w;

// %% eq(5) 
exp(X2ww_ss)*(X2ww - X2ww_ss) = (1/(exp(C_ss)*(1-chic/exp(dA_ss)))*(exp(W_ss)^epsw)*exp(N_ss))*((-1/(1-chic/exp(dA_ss))*((C - C_ss) - chic/exp(dA_ss)*(C(-1) - C_ss) + chic/exp(dA_ss)*(dA - dA_ss)) + epsw*(W - W_ss) + (N - N_ss))) 
+ bet*thew*(exp(dA_ss)^(epsw-1))*exp(X2ww_ss)*((epsw-1)*(dA(+1) - dA_ss) + (X2ww(+1) - X2ww_ss)) ; 

// %% eq(6)                                                                                                            
(w-w_ss) = (Zw-Zw_ss) + (X1w-X1w_ss) - (X2w-X2w_ss);

// %% eq(7) 식 수정 필요
exp(W)   = (thew*exp(W(-1)/exp(dA))^(1-epsw) + (1-thew)*exp(w)^(1-epsw))^(1/(1-epsw));
// exp(W)   = (thew*(exp(W(-1))/exp(dA))^(1-epsw) + (1-thew)*exp(w)^(1-epsw))^(1/(1-epsw));                                                    

/////////////////////////// 소비복합재 관련 조건식 ///////////////////////////////////////////

// %% eq(8)
(Cd - Cd_ss) = rho_c*(Pc - Pc_ss) + (C - C_ss) ;

// %% eq(9)
(Cf - Cf_ss) = rho_c*((Pc - Pc_ss) - (PPf - PPf_ss)) + (C - C_ss) ;

// %% eq(10) 
(1-rho_c)*((exp(Pc_ss))^(1-rho_c))*(Pc - Pc_ss) = (1-alp_c)*((exp(PPf_ss))^(1-rho_c))*(1-rho_c)*(PPf - PPf_ss)  ;

// %% eq(11)               
dP + Zinf = Pc - Pc(-1) + dPd ;

///////////////////////////// 가격 경직성 조건식 (국내 및 수입 소비재) ////////////////////////////////////////////////////

// %% eq(12) Log-linear 필요 (Zmc 충격관련 변수 처리방법)
//exp(X1)         = exp(Y)*exp(MC) + the*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *exp(dPd(+1))^(Zmc+1)*exp(dPd)^((iota)*(-Zmc)) *pi_target^((-Zmc)*(1-iota))*exp(X1(+1))  ;

 (exp(MC_ss)*exp(Y_ss) + the*bet*(exp(dPd_ss)^(Zmc + 1 - Zmc*iota)*pi_target^(-Zmc + Zmc*iota)*exp(X1_ss)/exp(dP_ss)))*(X1 - X1_ss) = exp(MC_ss)*exp(Y_ss)*((MC - MC_ss) + (Y - Y_ss))
 + the*bet*exp(dPd_ss)^(Zmc + 1 - Zmc*iota)*pi_target^(-Zmc + Zmc*iota)*exp(X1_ss)/exp(dP_ss)*
 ((1/(1 - chic/exp(dA_ss)))*((C - C_ss) - chic/exp(dA_ss)*(C(-1) - C_ss) + chic/exp(dA_ss)*(dA(-1) - dA_ss) - (C(+1) - C_ss) + chic/exp(dA_ss)*(C - C_ss) - chic/exp(dA_ss)*(dA - dA_ss))
 + (Zmc + 1)*(dPd(+1) - dPd_ss) - Zmc*iota*(dPd - dPd_ss) + (X1(+1) - X1_ss) - (dP(+1) - dP_ss)) ; 

// %% eq(13) Log-linear 필요 (Zmc 충격관련 변수 처리방법)
//exp(X2)         = exp(Y)         + the*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *exp(dPd(+1))^(Zmc)  *exp(dPd)^((iota)*(1-Zmc))*pi_target^((1-Zmc)*(1-iota))*exp(X2(+1))  ;
//

 (exp(Y_ss) + the*bet*(exp(dPd_ss)^(Zmc + iota - Zmc*iota)*pi_target^(1-iota-Zmc + Zmc*iota)*exp(X2_ss)/exp(dP_ss)))*(X2 - X2_ss) = exp(Y_ss)*(Y - Y_ss)
 + the*bet*exp(dPd_ss)^(Zmc + iota - Zmc*iota)*pi_target^(1-iota-Zmc + Zmc*iota)*exp(X2_ss)/exp(dP_ss)*
 ((1/(1 - chic/exp(dA_ss)))*((C - C_ss) - chic/exp(dA_ss)*(C(-1) - C_ss) + chic/exp(dA_ss)*(dA(-1) - dA_ss) - (C(+1) - C_ss) + chic/exp(dA_ss)*(C - C_ss) - chic/exp(dA_ss)*(dA - dA_ss))
 + Zmc*(dPd(+1) - dPd_ss) + (1-Zmc)*iota*(dPd - dPd_ss) + (X2(+1) - X2_ss) - (dP(+1) - dP_ss)) ; 



//exp(X2_ss)*(X2-X2_ss) = exp(Y_ss)*(Y-Y_ss)+(the*bet*(1/exp(dP_ss))*exp(dPd_ss)^(Zmc+iota*(1-Zmc))*pi_target^((1-Zmc)*(1-iota))*exp(X2_ss))
//* ( (1/(1-chic/exp(dA_ss)))*((C-C_ss)-((C(-1)-C_ss)-(dA(-1)-dA_ss))*(chic/exp(dA_ss))-(C(+1)-C_ss)+((C-C_ss)-(dA-dA_ss))*(chic/exp(dA_ss)))
//- (dP(+1)-dP_ss) + Zmc*(dPd(+1)-dPd_ss)+iota*(1-Zmc)*(dPd-dPd_ss)+(X2(+1)-X2_ss));

// %% eq(14)
(P - P_ss) = (X1 - X1_ss) - (X2 - X2_ss) ;

// %% eq(15) Log-linear 필요 (Zmc 충격관련 변수 처리방법)
1               = exp(dPd(-1))^(iota*(1-Zmc))*pi_target^((1-iota)*(1-Zmc))*the*exp(dPd)^(Zmc-1) + (1-the)*exp(P)^(1-Zmc);
// 0 = the*(dPd_ss)^((1-iota)*(Zmc-1))*pi_target^((1-iota)*(1-Zmc))*(iota*(1-Zmc)*dPd(-1) + (Zmc-1)*dPd(-1)) + ((1-the)*(P_ss^(1-Zmc))*(1-Zmc)*P) ;

// %% eq(16) Log-linear 필요 : Price Distortion
exp(SS)         = (exp(dPd(-1))^iota*pi_target^(1-iota))^(-eps)*the*exp(dPd)^(eps)*exp(SS(-1)) + (1-the)*exp(P)^(-eps);
//SS_ss*SS = the*SS_ss*SS(-1)*(pi_target^(-eps*(1-iota)))*(dPd_ss^(-eps*(iota-1)))*((-eps)*iota*dPd(-1) + eps*dPd) - (1-the)*(P_ss^(-eps))*eps*P ;

// %% eq(17) Log-linear 필요
//exp(X1f)         = exp(Cf)*exp(ERE)*exp(Pc)*(1+v*(exp(R)-1))*exp(Pistar) + thef*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *(exp(dPf(+1)))^(epsf)*(exp(dPf))^((iota)*(-epsf)) *pi_target^(-epsf*(1-iota))*exp(dPd(+1)) *exp(X1f(+1))  ;
// X1f_ss*X1f = (ERE_ss*Pc_ss*(1+v*(exp(R_ss)-1))*exp(Pistar_ss))*Cf_ss*(ERE_ss

(exp(Cf_ss)*exp(ERE_ss)*exp(Pc_ss)*(1+v*(exp(R_ss)-1))*exp(Pistar_ss) + thef*bet*(exp(dPf_ss))^(epsf-iota*epsf)*pi_target^(-epsf*(1-iota))*exp(dPd_ss)*exp(X1f_ss)/exp(dP_ss)) *(X1f - X1f_ss) =
exp(Cf_ss)*exp(ERE_ss)*exp(Pc_ss)*(1+v*(exp(R_ss)-1))*exp(Pistar_ss)*( (Cf-Cf_ss) + (ERE-ERE_ss) + (Pc-Pc_ss) + (Pistar - Pistar_ss) + (v*exp(R_ss)/(1+v*(exp(R_ss)-1))) * (R - R_ss))
+ thef*bet*(exp(dPf_ss))^(epsf-iota*epsf)*pi_target^(-epsf*(1-iota))*exp(dPd_ss)*exp(X1f_ss)/exp(dP_ss)*
((1/(1 - chic/exp(dA_ss)))*((C - C_ss) - chic/exp(dA_ss)*(C(-1) - C_ss) + chic/exp(dA_ss)*(dA(-1) - dA_ss) - (C(+1) - C_ss) + chic/exp(dA_ss)*(C - C_ss) - chic/exp(dA_ss)*(dA - dA_ss))
+ epsf*(dPf(+1)-dPf_ss) + iota*(-epsf)*(dPf-dPf_ss) +(dPd(+1)-dPd_ss)+(X1f(+1)-X1f_ss)-(dP(+1)-dP_ss));

// %% eq(18) Log-linear 필요
//exp(X2f)         = exp(Cf)                                               + thef*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *(exp(dPf(+1)))^(epsf)*(exp(dPf))^((iota)*(1-epsf)) *pi_target^((1-epsf)*(1-iota)) *exp(X2f(+1))  ;
//

(exp(Cf_ss) + thef*bet*exp(dPf_ss)^(epsf+iota*(1-epsf))*pi_target^((1-epsf)*(1-iota)) *exp(X2f_ss)/exp(dP_ss)) *(X2f - X2f_ss) = 
exp(Cf_ss)*(Cf-Cf_ss) + thef*bet*(exp(dPf_ss))^(epsf+iota*(1-epsf))*pi_target^((1-epsf)*(1-iota))*exp(X2f_ss)/exp(dP_ss)*
((1/(1 - chic/exp(dA_ss)))*((C - C_ss) - chic/exp(dA_ss)*(C(-1) - C_ss) + chic/exp(dA_ss)*(dA(-1) - dA_ss) - (C(+1) - C_ss) + chic/exp(dA_ss)*(C - C_ss) - chic/exp(dA_ss)*(dA - dA_ss))
+ epsf*(dPf(+1)-dPf_ss) + iota*(1-epsf)*(dPf-dPf_ss) + (X2f(+1) - X2f_ss) - (dP(+1) - dP_ss));

// %% eq(19)
(Pf - Pf_ss) = (X1f - X1f_ss) - (X2f - X2f_ss) ;

// %% eq(20) 
0 = thef*(exp(dPf_ss)^((1-iota)*(epsf-1)))*(pi_target^((1-iota)*(1-epsf)))*(iota*(1-epsf)*(dPf(-1)-dPf_ss) + (epsf-1)*(dPf-dPf_ss)) + (1-thef)*((exp(Pf_ss)/exp(PPf_ss))^(1-epsf))*(1-epsf)*((Pf-Pf_ss) - (PPf-PPf_ss)) ;

// %% eq(21)
 dPf + Zinff = PPf - PPf(-1) + dPd ;
 
 
/////////////////////////// 투자복합재 관련 조건식 ////////////////////////////////////////////

// %% eq(22)
(Invd - Invd_ss) = rho_i*(Pinv - Pinv_ss) + (Inv - Inv_ss) ;

// %% eq(23)
(Invf - Invf_ss) = rho_i*((Pinv - Pinv_ss) - (PP_invf - PP_invf_ss)) + (Inv - Inv_ss) ;

// %% eq(24) 
(1-rho_i)*((exp(Pinv_ss)^(1-rho_i))*(Pinv - Pinv_ss)) = (1-alp_i)*(exp(PP_invf_ss)^(1-rho_i))*(1-rho_i)*(PP_invf - PP_invf_ss)  ;

// %% eq(25)             
 dPinv = Pinv - Pinv(-1) + dPd ;

// %% eq(26) Log-linear 필요
exp(X1_invf)    = exp(Invf)*exp(ERE)*exp(Pc)*(1+v*(exp(R)-1))*exp(Pistar) + thef*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *(exp(dP_invf(+1)))^(epsf) *(exp(dP_invf))^((iota)*(-epsf)) *pi_target^(-epsf*(1-iota))*exp(dPd(+1))*exp(X1_invf(+1))  ;
//X1_invf_ss*(X1_invf - X1_invf_ss) = exp(Invf_ss)*(exp(ERE_ss)*exp(Pc_ss)*(1+v*(exp(R_ss)-1))*exp(Pistar_ss))*((ERE - ERE_ss) + (Pc - Pc_ss) + (v*(exp(R_ss)/(1+v*(exp(R_ss)-1))))*(R - R_ss) + (Invf - Invf_ss) + (Pistar - Pistar_ss))
//+ thef*bet*((exp(dP_invf_ss))^(epsf*(1-iota)))*(pi_target^(-epsf*(1-iota)))*exp(dPd_ss)*exp(X1_invf_ss)*((1/(1-chic/exp(dA_ss)))*(-(C(+1) - C_ss) + (chic/exp(dA_ss)+1)*(C - C_ss) - (chic/exp(dA_ss))*(C(-1) - C_ss) - (chic/exp(dA_ss))*(dA - dA_ss) + (chic/exp(dA_ss))*(dA(-1) - dA_ss)) + epsf*(dP_invf(+1) - dP_invf_ss) - epsf*iota*(dP_invf - dP_invf_ss) + (dPd(+1) - dPd_ss) + (X1_invf(+1) - X1_invf_ss)) ;

// %% eq(27) Log-linear 필요
exp(X2_invf)    = exp(Invf)                                               + thef*bet*((exp(C) - chic*exp(C(-1))/exp(dA(-1)) )/ (exp(C(+1)) - chic*exp(C)/exp(dA) ))/exp(dP(+1)) *(exp(dP_invf(+1)))^(epsf) *(exp(dP_invf))^((iota)*(1-epsf)) *pi_target^((1-epsf)*(1-iota))*exp(X2_invf(+1))  ;
//X2_invf - X2_invf_ss = (exp(Invf_ss)/exp(X2_invf_ss))*(Invf - Invf_ss) + thef*bet*(exp(dP_invf_ss))^(epsf+(1-epsf)*iota)*pi_target^((1-epsf)*(1-iota))*(
//(1/(1-chic/exp(dA_ss)))*(-(C(+1) - C_ss) + (chic/exp(dA_ss)+1)*(C - C_ss) - (chic/exp(dA_ss))*(C(-1) - C_ss) - (chic/exp(dA_ss))*(dA - dA_ss) + (chic/exp(dA_ss))*(dA(-1) - dA_ss))
//+ epsf*(dP_invf(+1)-dP_invf_ss) + (1-epsf)*iota*(dP_invf-dP_invf_ss)+(X2_invf(+1)-X2_invf_ss)) ;

// %% eq(28)
(P_invf - P_invf_ss) = (X1_invf - X1_invf_ss) - (X2_invf - X2_invf_ss) ;

// %% eq(29) 
 0 = thef*(exp(dP_invf_ss)^((1-iota)*(epsf-1))*pi_target^((1-iota)*(1-epsf)))*(iota*(1-epsf)*(dP_invf(-1) - dP_invf_ss) + (epsf-1)*(dP_invf - dP_invf_ss))
 + (1-thef)*(1-epsf)*((exp(P_invf)/exp(PP_invf))^(1-epsf))*((P_invf - P_invf_ss) - (PP_invf - PP_invf_ss)) ;

// %% eq(30) (좌변)의 Zinff를 곱하기에서 나누기로 변경
 dP_invf - Zinff = PP_invf - PP_invf(-1) + dPd ;


//////////////////////////////////// 금융가속기 (Financial Accelerator) ////////////////////////////////////////////

// %% eq(31)
 exp(Re_ss)*(Re - Re_ss) = (exp(Rk_ss)*exp(ut_ss)*exp(dPd_ss)/exp(Q_ss))*((Rk - Rk_ss) + (ut - ut_ss) + (dPd - dPd_ss) - (Q(-1) - Q_ss)) 
 - ((gam*(exp(ut_ss)-1) + gam2/2*(exp(ut_ss)-1)^2)*exp(dPd_ss)/exp(Q_ss))*(((gam*(exp(ut_ss)) + gam2*exp(ut_ss)*(exp(ut_ss)-1))/(gam*(exp(ut_ss)-1) + gam2/2*(exp(ut_ss)-1)^2))*(ut - ut_ss) - (Zinv - Zinv_ss) + (dPd - dPd_ss) - (Q(-1) - Q_ss))
 + (1-del)*exp(dPd_ss)*((Q - Q_ss) + (dPd - dPd_ss) - (Q(-1) - Q_ss)) ;

// %% eq(32) 
//exp(NN)          = exp(Zn)*gamma1*( exp(Re)*exp(Q(-1))*exp(K(-1))/exp(dA) - ( exp(NN(-1))/ (exp(K(-1))* exp(Q(-1))) )^(-omega)*exp(R(-1))/exp(dP)*(exp(Q(-1))*exp(K(-1))/exp(dA)- exp(NN(-1))/exp(dA)));
(NN - NN_ss) = (Zn - Zn_ss) + (1/(exp(NN_ss)/(exp(Zn_ss)*gamma1))) * (
  (exp(Re_ss)*exp(Q_ss)*exp(K_ss)/exp(dA_ss))*( (Re - Re_ss) + (Q(-1) - Q_ss) + (K(-1) - K_ss) - (dA - dA_ss))
+ omega*((exp(NN_ss)/(exp(Q_ss)*exp(K_ss)))^(-omega)) * (exp(R_ss)/exp(dP_ss)) * ((exp(Q_ss)*exp(K_ss)-exp(NN_ss))/exp(dA_ss)) * ( (NN(-1)-NN_ss) - (Q(-1)-Q_ss) - (K(-1)-K_ss))
- ((exp(NN_ss)/(exp(Q_ss)*exp(K_ss)))^(-omega)) * (exp(R_ss)/exp(dP_ss)) * ((exp(Q_ss)*exp(K_ss)-exp(NN_ss))/exp(dA_ss)) * ((R(-1)-R_ss) - (dP-dP_ss))
- ((exp(NN_ss)/(exp(Q_ss)*exp(K_ss)))^(-omega)) * (exp(R_ss)/exp(dP_ss)) * (
(exp(Q_ss)*exp(K_ss)/exp(dA_ss))*((Q(-1)-Q_ss)+(K(-1)-K_ss)) - (exp(NN_ss)/exp(dA_ss)) * (NN(-1)-NN_ss) - ((exp(Q_ss)*exp(K_ss)-exp(NN_ss))/exp(dA_ss)) * (dA - dA_ss)
)
) ;

// %% eq(33)
(Re(+1) - Re_ss) = (-omega)*((NN - NN_ss) - (Q - Q_ss) - (K - K_ss)) + (R - R_ss) - (dP(+1) - dP_ss) ;

// %% eq(34) 
 ( 1 - (1 - del)/exp(dA_ss))*((Zinv - Zinv_ss) + (Inv - Inv_ss)) + (exp(Zinv_ss)*exp(Inv_ss)/exp(K_ss))*(-chik*(exp(dA_ss) - exp(gamma)))*((Inv - Inv_ss) + (dA - dA_ss) - (Inv(-1) - Inv_ss))
= (K - K_ss) -  (1 - del)/exp(dA_ss)*((K(-1) - K_ss) - (dA - dA_ss)) ;    %% 투자 동태식 

// %% eq(35) Log_linear 새로고침
(exp(Rk_ss)*exp(Zinv_ss))*((Rk - Rk_ss) + (Zinv - Zinv_ss)) = (gam2*exp(ut_ss))*(ut - ut_ss) ;

// %% eq(36) 
 0 = exp(Zinv_ss)*exp(Q_ss)*(
 (1 - chik/2*(exp(dA_ss) - exp(gamma))^2 - chik*(exp(dA_ss) - exp(gamma))*dA_ss)*((Q - Q_ss) + (Zinv - Zinv_ss))
 +(-2*chik*(exp(dA_ss) - exp(gamma)) - chik*exp(dA_ss)^2)*((Inv - Inv_ss) + (dA - dA_ss) - (Inv(-1) - Inv_ss))
 )
 + bet*exp(Q_ss)*(
 chik*(exp(dA_ss) - exp(gamma))*dA_ss*((Q(+1) - Q_ss) - (dA(+1) - dA_ss))
 + (1/(1 - chic/exp(dA_ss)))*chik*(exp(dA_ss) - exp(gamma))*exp(dA_ss)*((C - C_ss) - chic/exp(dA_ss)*(C(-1) - C_ss) + chic/(exp(dA_ss)^2)*(dA - dA_ss) - (C(+1) - C_ss) + chic/exp(dA_ss)*(C - C_ss) - chic/(exp(dA_ss)^2)*(dA - dA_ss))
 + chik*(exp(dA_ss)^2)*((Inv(+1) - Inv_ss) + (dA(+1) - dA_ss) - (Inv - Inv_ss))
 + 2*bet*exp(Q_ss)*chik*(exp(dA_ss) - exp(gamma))*exp(dA_ss)*((Inv(+1) - Inv_ss) + (dA(+1) - dA_ss) - (Inv - Inv_ss))
 ) ;  %% 투자에 대한 일계조건식 (자본재 가격 결정식)



////////////////////////////////// 생산부문 //////////////////////////////////////////////////

// %% eq(37)                                                     
 (W - W_ss) - (Rk - Rk_ss) + exp(R_ss)/(1 + v*(exp(R_ss)-1))*(R - R_ss) = (K(-1) - K_ss) + (ut - ut_ss) - (N - N_ss) - (dA - dA_ss)  ;    // 자본 노동의 생산투입요소 비율

// %% eq(38) 
(MC - MC_ss) = -(Za - Za_ss) + (1-alpk)*(W - W_ss) + alpk*(Rk - Rk_ss) + (1-alpk)*v*exp(R_ss)*R/(1 + v*(exp(R_ss) - 1)) ;                // 마크업 정의식 

// %% 본 모형에서 제외된 식 
//exp(W)*(1+v*(exp(R)-1)) = exp(Za)*exp(MC)* (1-alpk)*(exp(K(-1))*exp(ut))^(alpk)   *exp(N)^(-alpk) *exp(dA)^(-alpk)/exp(SS);
//exp(Rk)                 = exp(Za)*exp(MC)* alpk*(exp(K(-1))*exp(ut))^(alpk-1)  *exp(N)^(1-alpk)*exp(dA)^(1-alpk)/exp(SS);

// %% eq(39) 
Y - Y_ss = (Za - Za_ss) + alpk*((K(-1) - K_ss) + (ut - ut_ss)) + (1-alpk)*(N - N_ss) - alpk*(dA - dA_ss) - (SS - SS_ss) ;

// %% eq(40) Log-linear 필요
exp(Y)                  = exp(Cd) + exp(Invd) + exp(Xd) + chik/2*(exp(Inv)*exp(dA)/exp(Inv(-1))-exp(gamma))^2*exp(Inv) + exp(G) + 1/exp(Zinv)*(gam*(exp(ut)-1) + gam2*(exp(ut)-1)^2)*exp(ut)*exp(K(-1));
//Y_ss*Y = Cd_ss*Cd + Invd_ss*Invd + Xd_ss*Xd + (dA_ss/(dA_ss - exp(gamma)))*(Inv + dA - Inv(-1)) + G_ss*G + (exp(ut_ss)*(exp(ut_ss)-1)*(gam + gam2*(exp(ut_ss)-1)/2))*ut + K_ss*K - Zinv_ss*Zinv  ;

// %% eq(41)  
(G - G_ss) = (Zg - Zg_ss) + (Y - Y_ss) - phi_g*(Y - steady_state(Y)) ;


/////////////////////////////////// 경제충격 프로세스 /////////////////////////////////////////////////////

// %% eq(42)
dA         = rhogam*dA(-1) + (1-rhogam)*gamma + etrend;                                                // 영구적 기술충격 ln(A_t) = ln(A_t-1) + gamma*exp(Zg) + ea //

// %% eq(43)
Za         = rhoa    *Za(-1)      + ea;

// %% eq(44)
Zinv       = rhoinv  *Zinv(-1)    + einv;

// %% eq(45)
Zc         = rhoc    *Zc(-1)      + ec;

// %% eq(46)
Zinf       = rhoinf  *Zinf(-1)    + einf;

// %% eq(47)
Zn         = rhon    *Zn(-1)      + en;

// %% eq(48)
Zinff      = rhoinff *Zinff(-1)   + einff;

// %% eq(49)
Zx         = rhox    *Zx(-1)      + ex;

// %% eq(50)
Zg         = rhog    *Zg(-1)      + eg;

// %% eq(51)
Zw         = rhow    *Zw(-1)      + ew;

// %% eq(52)
Zcp        = rhocp   *Zcp(-1)     + ecp;

// %% eq(53)
Zmc        = (1-rhomc)*eps + rhomc   *Zmc(-1)   + emc;


////////////////////////////////////// 통화준칙 (Taylor rule) /////////////////////////////////////////////////////////

// %% eq(54)
log(exp(R))     =   rhor*log(exp(R(-1)))+ (1-rhor)*(log(steady_state(exp(R)))+ (1+rhop)*log(exp(dP))+ rhoy*log(exp(Y)/steady_state(exp(Y)))) +er ;


////////////////////////////////////// Foreign VAR/////////////////////////////////////////////////////////

// %% eq(55)
Ystar           = a11*Ystar(-1)  + a12*Rstar(-1) + a13*Pistar(-1) + eystar;

// %% eq(56)
Rstar           = a21*Ystar(-1)  + a22*Rstar(-1) + a23*Pistar(-1) + erstar;

// %% eq(57)
Pistar          = a31*Ystar(-1)  + a32*Rstar(-1) + a33*Pistar(-1) + epistar;

 
//////////////////////////////////////// Export function & BOP  //////////////////////////////////////////

// %% eq(58)  
(X - X_ss) = (Zx - Zx_ss) + psi_x*(X(-1) - X_ss) - psi_x*(dA - dA_ss) - psi_px*(1 - psi_x)*((Pxx - Pxx_ss) - (S - S_ss) + (v*exp(R_ss)/(1 + v*(exp(R_ss) - 1)))*(R - R_ss)) + (1 - psi_x)*(Ystar - Ystar_ss) ;  // 수출함수

// %% eq(59)
Xd-Xd_ss = etax*(Pxx-Pxx_ss) + (X-X_ss) ; 

// %% eq(60) 
(Xm-Xm_ss) = etax*((Pxx-Pxx_ss) - (ERE-ERE_ss) - (Pc-Pc_ss)) + (X-X_ss) ;

// %% eq(61) 
(1 - etax)*(Pxx - Pxx_ss) = (1 - alpx)*(1 - etax)*((exp(ERE)*exp(Pc))^(1 - etax))/(alpx + (1 - alpx)*((exp(ERE)*exp(Pc))^(1 - etax)))*((ERE - ERE_ss) + (Pc - Pc_ss)) ;

// %% eq(62)
dPx = Pxx - Pxx(-1) + dPd ;

// %% eq(63) Log-linear 필요
exp(ERE)*exp(Px)*exp(Pc)*exp(X) + exp(S)*exp(RR(-1))*exp(F(-1))*exp(Pc)/(exp(dA)*exp(dP))  =  exp(S)*exp(F)*exp(Pc) + exp(ERE)*exp(Pc)*exp(Cf)*(1+v*(exp(R)-1)) + exp(ERE)*exp(Pc)*exp(Invf)*(1+v*(exp(R)-1)) + exp(ERE)*exp(Pc)*exp(Xm)*(1+v*(exp(R)-1));

//(exp(S_ss)*exp(F_ss)*exp(Pc_ss) + exp(ERE_ss)*exp(Pc_ss)*(1+v*(exp(R_ss)-1))*(exp(Cf_ss) + exp(Invf_ss) + exp(Xm_ss)))
//*(
//exp(ERE_ss)*exp(Px_ss)*exp(Pc_ss)*exp(X_ss)*((ERE - ERE_ss) + (Px - Px_ss) + (Pc - Pc_ss) + (X - X_ss))
//+ (exp(S_ss)*exp(RR_ss)*exp(F_ss)*exp(Pc_ss)/(exp(dA_ss)*exp(dP_ss))*((S - S_ss) + (RR(-1) - RR_ss) + (F(-1) - F_ss) + (Pc - Pc_ss) - (dA - dA_ss) - (dP - dP_ss)))
//) 
//= (exp(ERE_ss)*exp(Px_ss)*exp(Pc_ss)*exp(X_ss) + exp(S_ss)*exp(RR_ss)*exp(F_ss)*exp(Pc_ss)/(exp(dA_ss)*exp(dP_ss)))
//*(
//(exp(S_ss)*exp(F_ss)*exp(Pc_ss))*((S - S_ss) + (F - F_ss) + (Pc - Pc_ss))
//+ exp(ERE_ss)*exp(Pc_ss)*(exp(Cf_ss) + exp(Invf_ss) + exp(Xm_ss))*((1+v*(exp(R_ss)-1))*((ERE - ERE_ss) + (Pc - Pc_ss)) + v*exp(R_ss)*(R - R_ss))
//+ exp(ERE_ss)*exp(Pc_ss)*(1+v*(exp(R_ss)-1))*(exp(Cf_ss)*(Cf - Cf_ss) + exp(Invf_ss)*(Invf - Invf_ss) + exp(Xm_ss)*(Xm - Xm_ss))
//) ;




////////////////////////////////////////  UIP condition //////////////////////////////////////////

// %% eq(64)
(R - R_ss) = (RR - RR_ss) + (S(+1) - S) ;                                                     // UIP condition

// %% eq(65)
//exp(RR)         = exp(RRw)*(exp(Rstar)) * exp(-zet*(exp(S)*exp(F)- (exp(steady_state(Y))*phi)) - zet2*(exp(RR)-exp(steady_state(RR))- (exp(R)-exp(steady_state(R)))))*exp(Zcp);                     // 

(RR - RR_ss) = (RRw - RRw) + (Rstar - Rstar_ss) - zet*(exp(S_ss)*exp(F_ss))*((S - S_ss) + (F - F_ss))  - zet2*exp(RR_ss)*(RR - RR_ss) + zet2*exp(R_ss)*(R - R_ss) + (Zcp - Zcp_ss) ;

//exp(RR_ss)*(RR - RR_ss) = (exp(RRw_ss)*exp(Rstar_ss)*exp(-zet*(exp(S_ss)*exp(F_ss)- (exp(Y_ss)*phi)))*exp(Zcp_ss))
//*((RRw - RRw_ss) + (Rstar - Rstar_ss) + (Zcp - Zcp_ss) 
//+ (-zet*(exp(S)*exp(F) -exp(Y_ss)*phi)-zet2*(exp(RR)-exp(RR_ss)-(exp(R)-exp(R_ss))))
//- (-zet*(exp(S_ss)*exp(F_ss)-exp(Y_ss)*phi)));
 
 
 
// %% eq(66)
RRw = RRw_ss ;

// %% eq(67)
Px = PPx ;

// %% eq(68)
Pm = PPm + Pistar - S ;                                                                    // 수입재 가격

// %% eq(69)
ERE = PPm + Pistar + S - Pc;                                                              // 실질환율

// %% eq(70)
PPx = PPx_ss;

// %% eq(71)
PPm = PPm_ss;

// %% eq(72)
Imp*Imp_ss = Cf*Cf_ss + Invf*Invf_ss + Xm*Xm_ss ;
  


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

//observation_trends;
//Y_obs (gamma);
//C_obs (gamma);
//Inv_obs (gamma);
//X_obs (gamma);
//G_obs (gamma);
//W_obs (gamma);

//end;


initval;
Y               	=	1.13119084748	;
C               	=	0.60625951214	;
N               	=	0.37369265825	;
R               	=	0.01050050033	;
W               	=	0.17811217949	;
Rk              	=	-2.81721866732	;
K               	=	2.64823304686	;
dA              	=	0.00950000000	;
dP              	=	0.00000000000	;
G               	=	-0.47824706496	;
MC              	=	-0.16913302902	;
P               	=	-0.01361141238	;
X1              	=	2.51625764962	;
X2              	=	2.69821987121	;
SS              	=	0.00072508622	;
ERE             	=	0.04766369226	;
F               	=	-0.84031646611	;
RR              	=	0.01050050033	;
S               	=	0.16690158463	;
Px              	=	0.00000172120	;
Pm              	=	0.16690198542	;
X               	=	0.01746058156	;
PPx             	=	0.00000040050	;
PPm             	=	0.00000040050	;
RRw             	=	0.01045766585	;
Q               	=	0.00000000002	;
Inv             	=	-0.72675164538	;
Re              	=	0.03481780798	;
NN              	=	2.16188689389	;
Cd              	=	0.51085654666	;
Cf              	=	-1.31090409273	;
Pc              	=	0.11923829287	;
Pf              	=	0.38817464729	;
dPd             	=	0.00000000000	;
dPf             	=	0.00000000008	;
PPf             	=	0.39221294140	;
PP_invf         	=	0.39221240374	;
X1f             	=	-0.36582136226	;
X2f             	=	-0.53927102871	;
Invd            	=	-1.36832585679	;
Invf            	=	-1.40932389346	;
Pinv            	=	0.24252049187	;
P_invf          	=	0.38817446823	;
dPinv           	=	0.00000000000	;
dP_invf         	=	0.00000044990	;
X1_invf         	=	-0.46423932484	;
X2_invf         	=	-0.63768881222	;
Imp             	=	0.00744009976	;
Ystar           	=	0.00000000000	;
Rstar           	=	0.00000000000	;
Pistar          	=	0.00000000000	;
Za              	=	0.00000000000	;
Zinv            	=	0.00000000000	;
Zc              	=	0.00000000000	;
Zinf            	=	0.00000000000	;
Zn              	=	0.00000000000	;
Zinff           	=	0.00000000000	;
Zx              	=	0.00000000000	;
Zg              	=	0.00000000000	;
Zw              	=	0.00000000000	;
Zcp             	=	0.00000000000	;
Zmc             	=	6.45400000000	;
Y_obs           	=	0.00000000000	;
R_obs           	=	0.00000000000	;
dP_obs          	=	0.00000000000	;
C_obs           	=	0.00000000000	;
Inv_obs         	=	0.00000000000	;
NN_obs          	=	0.00000000000	;
dPf_obs         	=	0.94000000822	;
X_obs           	=	0.00000000000	;
dPinv_obs       	=	0.00000000000	;
S_obs           	=	0.00000000000	;
Ystar_obs       	=	0.00000000000	;
Rstar_obs       	=	0.00000000000	;
Pistar_obs      	=	0.00000000000	;
W_obs           	=	0.00000000000	;
G_obs           	=	0.00000000000	;
SP_obs          	=	0.00000000000	;
dPx_obs         	=	-0.25800000000	;
GAPeff          	=	0.00000000000	;
Y_e             	=	-0.02722764861	;
Ce              	=	-0.40167832198	;
Ne              	=	-1.10965997753	;
LAMe            	=	0.40167834884	;
We              	=	0.69145476234	;
R_e             	=	0.01050251301	;
Ke              	=	2.17993798893	;
Inve            	=	-1.19504670042	;
Rke             	=	-3.30632826186	;
Qe              	=	0.03031681247	;
w               	=	0.18011580569	;
X1w             	=	1.74675667913	;
X1ww            	=	4.16034449855	;
X2w             	=	1.74597887036	;
X2ww            	=	2.84306424283	;
ut              	=	0.03260590061	;
Xd              	=	-0.63618613273	;
Xm              	=	-0.70606406092	;
Pxx             	=	0.08918874193	;
dPx             	=	0.00000000000	;
Ystar_or        	=	0.00000000000	;

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

//varobs Y_obs, R_obs, dP_obs, C_obs, Inv_obs, NN_obs, X_obs,S_obs, W_obs, Ystar_obs, Rstar_obs, Pistar_obs, G_obs, dPx_obs;//, dPf_obs;//,  G_obs;//;// dPf_obs;//, W_obs ;//dPf_obs,  dPinv_obs ;//,  Ystar_obs, Rstar_obs, Pistar_obs G_obs, W_obs, dPd_obs;//, N_obs;//, SP_obs;

//estimated_params;
//    eta,        gamma_pdf,      1.4,     0.2;
//    eps,        gamma_pdf,      6,       1.0;
//    epsf,       gamma_pdf,      6,       1.0;
//    epsw,       gamma_pdf,      6,       1.0;
//    the,        beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05 
//    thef,       beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05
//    thew,       beta_pdf,       0.6,    0.05;   //0.70,  0.05;  % previously 0.7, 0.05
//    chic,       normal_pdf,     0.5,     0.1;
//    chik,       normal_pdf,     4,       1.0;
//    omega,      beta_pdf,       0.3,     0.05;    
    
//    rhor,       beta_pdf,       0.8,     0.05;
//    rhoy,       beta_pdf,       0.2,     0.05;
//    rhop,       beta_pdf,       0.5,     0.1;   
//    iota,       beta_pdf,       0.5,     0.1;
//    psi_x,      beta_pdf,       0.5,     0.1;
//    psi_px,     gamma_pdf,      0.5,     0.1;
//    alp_c,      beta_pdf,       0.6,     0.1;
//    alp_i,      beta_pdf,       0.6,     0.1;     
//    rho_c,      gamma_pdf,      2,       0.5;
//    rho_i,      gamma_pdf,      2,       0.5;

//    alpx,       beta_pdf,       0.5,     0.1;
//    etax,       gamma_pdf,      2.5,     0.5;
//    gam2,       beta_pdf,       0.5,     0.1;     
//    zet2,       beta_pdf,       0.5,     0.1;     
//    phi_g,      gamma_pdf,      0.5,     0.1;
    
//    rhoa,       beta_pdf,       0.8,    0.1;
//    rhoc,       beta_pdf,       0.8,    0.1;
//    rhogam,     beta_pdf,       0.8,    0.1;
//    rhoinf,     beta_pdf,       0.8,    0.1;
//    rhoinv,     beta_pdf,       0.8,    0.1;
//    rhon,       beta_pdf,       0.8,    0.1;
//    rhoinff,    beta_pdf,       0.8,    0.1;
//    rhox,       beta_pdf,       0.8,    0.1;
//    rhopinv,    beta_pdf,       0.8,    0.1;
//    rhog,       beta_pdf,       0.8,    0.1;
//    rhow,       beta_pdf,       0.8,    0.1;
//    rhocp,      beta_pdf,       0.8,    0.1;

//    a11,        normal_pdf,     0.9,    0.1;
//    a12,        normal_pdf,     0.0,    0.1;
//    a13,        normal_pdf,     0.0,    0.1;
//    a21,        normal_pdf,     0.0,    0.1;
//    a22,        normal_pdf,     0.9,    0.1;
//    a23,        normal_pdf,     0.0,    0.1;
//    a31,        normal_pdf,     0.0,    0.1;
//    a32,        normal_pdf,     0.0,    0.1;
//    a33,        normal_pdf,     0.6,    0.1;

//    stderr er,      inv_gamma_pdf, 0.01,  2;
//    stderr ea,      inv_gamma_pdf, 0.01,  2;
//    stderr ec,      inv_gamma_pdf, 0.01,  2;
//    stderr etrend,  inv_gamma_pdf, 0.01,  2;
//    stderr einf,    inv_gamma_pdf, 0.01,  2;
//    stderr einv,    inv_gamma_pdf, 0.01,  2;
//    stderr en,      inv_gamma_pdf, 0.01,  2;
//    stderr einff,   inv_gamma_pdf, 0.01,  2;
//    stderr ex,      inv_gamma_pdf, 0.01,  2;
//    stderr eg,      inv_gamma_pdf, 0.01,  2;
//    stderr emc,     inv_gamma_pdf, 0.01,  2;
//    stderr eystar,  inv_gamma_pdf, 0.01,  2;
//    stderr erstar,  inv_gamma_pdf, 0.01,  2;
//    stderr epistar, inv_gamma_pdf, 0.01,  2;
//    stderr ew,      inv_gamma_pdf, 0.01,  2;
//    stderr ecp,     inv_gamma_pdf, 0.01,  2;
//    stderr epx,    inv_gamma_pdf, 0.01,  2;



//end;

//estimation(
//datafile     = data_201414, //data_201344_sa_soe_dmean.xls, //,
//first_obs    = 1, 
//mh_replic    = 50000,
//mode_file    = BOKDSGE2014_0623_TR_UIP_W_markup_mode,
//mh_nblocks   = 1,
//forecast     = 12, 
//nobs       = [5:50],
//mh_drop      = 0.4, 
//mh_jscale    = 0.20, 
//mode_compute = 4, 
//plot_priors  = 0,
//smoother, 
//endogenous_prior,
//,mode_check
//filtered_vars,
//nograph ) Y_obs, R_obs,  dP_obs, C_obs, Inv_obs, NN_obs, dPf_obs, X_obs, dPinv_obs, S_obs, Ystar_obs, Rstar_obs, Pistar_obs dA G_obs W_obs Za,  Zinv, Zc,  Zinf, Zn, Zinff, Zx, Zg, Zw, Zcp;

//shock_decomposition   Y_obs, C_obs, R_obs, dP_obs, Inv_obs, GAPeff;
     

            
     