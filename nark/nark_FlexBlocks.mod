/*  Euler_f: C_f, R_f, dP_f
exp(dZtr(+1))*exp(Zc)/(exp(C_f)-CHIc*exp(C_f(-1))/exp(dZtr))
= BETA*exp(Zc(+1))/(exp(C(+1))-CHIc*exp(C)/exp(dZtr(+1)))*exp(R_f)/exp(dP_f(+1));
*/

/*  W_f : W_f, C_f
exp(W_f) = PSIw/(PSIw-1)*exp(Zw)*(exp(C_f)-CHIc*exp(C_f(-1))/exp(dZtr))*exp(N_f)^ETA;
*/

/*  Consumption Composite_f : Ch_f, P_f, C_f, Cf_f, Pcf_f, dP_f, dPh_f
exp(Ch_f) = ALPHAc*exp(P_f)^XIc*exp(C_f);
exp(Cf_f) = ((1-ALPHAc)*(exp(P_f)/exp(Pcf_f))^XIc*exp(C_f));
exp(P_f) = (ALPHAc + (1-ALPHAc)*(exp(Pcf_f))^(1-XIc))^(1/(1-XIc));
exp(dP_f)*exp(Zpi) = (exp(P_f)/exp(P_f(-1))*exp(dPh_f));
*/

/*  Ph_f    : Ph_f, MC_f, dPh_f
exp(Ph_f) = exp(MC_f);
exp(dPh_f) = exp(Ph_f)/exp(Ph_f(-1));
*/

/*Pcf_f   : Pcf_f, dPcf_f, q_f, P_f, R_f
exp(Pcf_f) = exp(q_f)*exp(P_f)*(1+v*(exp(R_f)-1));
exp(dPcf_f) = exp(Pcf_f)/exp(Pcf_f(-1));
*/

/*  I_f and Pi_f     : Ih_f, If_f, Pi_f, dPi_f, I_f, Pif_f, dPi_f, 
exp(Ih_f) = ALPHAi*exp(Pi_f)^XIi*exp(I_f);
exp(If_f) = (1-ALPHAi)*(exp(Pi_f)/exp(Pif_f))^XIi*exp(I_f);
exp(Pi_f) = (ALPHAi + (1-ALPHAi)*(exp(Pif_f))^(1-XIi))^(1/(1-XIi));
exp(dPi_f) = exp(Pi_f)/exp(Pi_f(-1))*exp(dPh_f);

exp(Pif_f) = exp(q_f)*exp(P_f)*(1+v*(exp(R_f)-1))
exp(dPif_f) = exp(Pif_f)/exp(Pif_f(-1));
*/

/* Financial Accelerator_f: Re_f, Q_f, Rk_f, ut_f, dPh_f, NW_f, K_f, R_f, dP_f, I_f
exp(Re_f) = (exp(Q_f)*(1-DELTA) + exp(Rk_f)*exp(ut_f)
             - 1/exp(Zi)*(GAMMA1*(exp(ut_f)-1)+GAMMA2/2*(exp(ut_f)-1)^2))
            *exp(dPh_f)/exp(Q_f(-1));
exp(NW_f) = exp(Znw)*SP
            *(exp(Re_f)*exp(Q_f(-1))*exp(K_f(-1))/exp(dZtr)
              - (exp(NW_f(-1))/(exp(K_f(-1))*exp(Q_f(-1))))^(-KAPPArp)*exp(R_f(-1))
                /exp(dP_f)*(exp(Q_f(-1))*exp(K_f(-1))/exp(dZtr)
                            - exp(NW_f(-1))/exp(dZtr)));
exp(Re_f(+1)) = (exp(NW_f)/(exp(K_f)*exp(Q_f)))^(-KAPPArp)*exp(R_f)/exp(dP_f(+1));

exp(K_f)-(1-DELTA)*exp(K_f(-1))/exp(dZtr)
= exp(Zi)*exp(I_f)*(1-KAPPAi/2*(exp(I_f)*exp(dZtr)/exp(I_f(-1))-exp(GAMMAtr))^2);

exp(Rk_f) = (GAMMA1+GAMMA2*(exp(ut_f)-1))/exp(Zi);

1 = exp(Q_f)*exp(Zi)*(1 - KAPPAi/2*(exp(I_f)*exp(dZtr)/exp(I_f(-1)) - exp(GAMMAtr))^2
                      - KAPPAi*(exp(I_f)*exp(dZtr)/exp(I_f(-1))-exp(GAMMAtr))
                        *exp(I_f)*exp(dZtr)/exp(I_f(-1)))
    + BETA*exp(Q_f(+1))*exp(Zi(+1))
      *(1/exp(dZtr(+1))*(exp(C_f)-CHIc*exp(C_f(-1))/exp(dZtr))
        /(exp(C_f(+1))-CHIc*exp(C_f)/exp(dZtr(+1))))
      *KAPPAi*(exp(I_f(+1))*exp(dZtr(+1))/exp(I_f)-exp(GAMMAtr))
      *(exp(I_f(+1))*exp(dZtr(+1))/exp(I_f))^2;
*/

/*  Production: W_f, R_f, Rk_f, K_f, ut_f, N_f, MC_f, Y_f
// MPK = MPL
exp(W_f)*(1+v*(exp(R_f)-1))/exp(Rk_f)
= (1-ALPHAk)/ALPHAk*exp(K_f(-1))*exp(ut_f)/exp(N_f)/exp(dZtr);

// MC
exp(MC_f) = 1/exp(Za)*(exp(W_f)*(1+v*(exp(R_f)-1)))^(1-ALPHAk)*exp(Rk_f)^ALPHAk
            *(1-ALPHAk)^(ALPHAk-1)*ALPHAk^(-ALPHAk);

exp(Y_f) = exp(Za)*(exp(K_f(-1))*exp(ut_f))^ALPHAk
           *exp(N_f)^(1-ALPHAk)*exp(dZtr)^(-ALPHAk);
exp(Y_f) = exp(Ch_f) + exp(Ih_f) + exp(Xh_f) + exp(G_f)
         + KAPPAi/2*(exp(I_f)*exp(dZtr)/exp(I_f(-1))-exp(GAMMAtr))^2*exp(I_f)
         + 1/exp(Zi)*(GAMMA1*(exp(ut_f)-1)+GAMMA2*(exp(ut_f)-1)^2)
           *exp(ut_f)*exp(K_f(-1));

exp(G_f) = Gss*exp((Y_f))*exp(Zg)*(exp(Y_f)/exp(steady_state(Y_f)))^(-PHIg);
*/

/*  Export function & BOP   :X_f, R_f, Px_f, S_f, Xh_f, Xf_f, q_f, P_f, dPx_f, dPh_f
exp(X_f) = exp(Zx)*(exp(X_f(-1))/exp(dZtr))^PHIx
           *(((1+v*(exp(R_f)-1))*exp(Px_f)/exp(S_f)/exp(PIstar))^(-XIh)
             *(exp(Ystar)))^(1-PHIx);
exp(Xh_f) = ALPHAx*exp(Px_f)^XIx*exp(X_f);
exp(Xf_f) = (1-ALPHAx)*(exp(Px_f)/(exp(q_f)*exp(P_f)))^XIx*exp(X_f);
exp(Px_f) = (ALPHAx+(1-ALPHAx)*(exp(q_f)*exp(P_f))^(1-XIx))^(1/(1-XIx));
exp(dPx_f)*exp(epx) = exp(Px_f)/exp(Px_f(-1))*exp(dPh_f);

exp(Px_f)*exp(X_f)
+ exp(S_f)*exp(Rf_f(-1))*exp(F_f(-1))*exp(P_f)/exp(dZtr)/exp(dP_f)
= exp(S_f)*exp(F_f)*exp(P_f)
  + exp(q_f)*exp(P_f)*(1+v*(exp(R_f)-1))*(exp(Cf_f) + exp(If_f) + exp(Xf_f));
*/

/*  UIP condition   : R_f, Rf_f, S_f, F_f, Y_f, q_f, M_f
exp(R_f) = exp(Rf_f)*exp(S_f(+1))/(exp(S_f));
exp(Rf_f) = exp(- UIPf*(exp(S_f)*exp(F_f)-(exp(steady_state(Y_f))*UIPy))
                - UIPr*(exp(Rf_f)-exp(steady_state(Rf_f))
                        - (exp(R_f)-exp(steady_state(R_f)))))
          *exp(Zcp)*exp(Rw)*(exp(Rstar));

exp(Rw) = (1+GAMMAtr)/BETA;                   
//exp(q) = exp(PPm)*(exp(PIstar))*exp(S)/exp(P);                          //??
//exp(q) = exp(q(-1))*(exp(PIstar))*exp(S)/exp(S(-1))*exp(P(-1))/exp(P);    //
//exp(q) = exp(Rf)*exp(q(+1))*exp(dP(+1))/exp(PIstar(+1))/exp(R);           //
exp(q_f) = exp(S_f)*exp(PIstar)/exp(P_f);
exp(M_f) = exp(Cf_f) + exp(If_f) + exp(Xf_f);
*/