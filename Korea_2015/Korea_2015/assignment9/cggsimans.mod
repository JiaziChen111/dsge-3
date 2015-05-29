// Simulation of the Clarida, Gali, Gertler model

var da pie r rstar tau x dy; 	   // DECLARATION OF THE ENDOGENOUS VARIABLES 

varexo  eps_a eps_tau ;    // DECLARATION OF THE STRUCTURAL INNOVATIONS.

parameters phi beta kappa phi_x phi_pie alpha rho lambda; // DECLARATION OF THE DEEP PARAMETERS.

// Parameter Values  
beta    = 0.97;
phi_x   = .15;
phi_pie = 1.5;
alpha   = 0.8;
rho     = 0.9;
lambda  = 0.5;
phi     = 1;
theta   = 0.75;
kappa   = ((1-theta)*(1-beta*theta)*(1+phi))/(theta);

// DECLARATION OF THE (LINEAR) DSGE MODEL: 
model(linear);

    beta*pie(+1) + kappa*x = pie;                                      // Calvo Pricing Equation
    r - pie(+1)-rstar   = x(+1) - x;                                   // Intertemporal Equation
    alpha*r(-1) + (1-alpha)*phi_pie*pie + (1-alpha)*phi_x*x = r;   // Monetary Policy Rule
    rstar = da(+1) + (1-lambda)/(1+phi)*tau;                           // Definition of the Natural Rate
    da = rho*da(-1) + eps_a;
    tau = lambda*tau(-1) + eps_tau;
    dy  = x - x(-1) + da - (tau - tau(-1))/(1+phi);                    // construction of output growth

end;


shocks;
var eps_a;
stderr .02;
var eps_tau;
stderr .02;
end;

set_dynare_seed=1;
stoch_simul(periods=5000, irf=7, nograph);

// answer to 2a:
//plots;

// partial answer to 2f
//fprintf('rstar = %5.6f, x = %5.6f, r = %5.6f, pie = %5.6f\n',rstar_signal(1), ...
//  x_signal(1), r_signal(1), pie_signal(1));


% answer to 3b:
figure
dysim=oo_.endo_simul(7,:);
gap=oo_.endo_simul(6,:);
ysim=cumsum(oo_.endo_simul(7,:));
lam=1;
[y_hp,y_hptrend]=HPFAST(ysim,lam);
tt=1:length(ysim);
tt=220:300;
subplot(221)
plot(tt,y_hptrend(tt),tt,ysim(tt))
title('lam=1')
axis tight;
lam=1600;
[y_hp,y_hptrend]=HPFAST(ysim,lam);
subplot(222);
plot(tt,y_hptrend(tt),tt,ysim(tt));
title('lam=1600');
axis tight;
lam=50000;
[y_hp,y_hptrend]=HPFAST(ysim,lam);
subplot(223);
plot(tt,y_hptrend(tt),tt,ysim(tt));
title('lam=50000');
axis tight;
lam=160000000;
[y_hp,y_hptrend]=HPFAST(ysim,lam);
subplot(224)
plot(tt,y_hptrend(tt),tt,ysim(tt))
title('lam=160000000')
axis tight;

figure
%answer to 3b:
subplot(211)
lam=1600;
%lam=50;
[y_hp,y_hptrend]=HPFAST(ysim,lam);
ystar=ysim-gap;

% computing the correlation, as suggested in the question...
ca=corrcoef(gap,y_hp);

% computing standard deviations...
sgap=std(gap);
syhp=std(y_hp);
plot(tt,y_hptrend(tt),tt,ystar(tt),'+-',tt,ysim(tt));
legend('hp trend','natural output','actual output')
axis tight;
subplot(212)
plot(tt,y_hp(tt),tt,gap(tt),'+-')
legend('hp-filtered output','actual gap')

% note how the correlations and standard deviations (i.e., numbers) are inserted into a string for printing:
str=['correlation, hp-filtered output and actual output gap = ',num2str(ca(1,2)),' std(gap) = ',num2str(sgap),' std(yhp) = ',num2str(syhp)];
title(str)
axis tight;

sim = oo_.endo_simul;

save data sim;