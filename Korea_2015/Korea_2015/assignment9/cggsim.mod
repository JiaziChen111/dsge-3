// Simulation of the Clarida, Gali, Gertler model

var a pie r rstar tau x; 	   // DECLARATION OF THE ENDOGENOUS VARIABLES

varexo eps_a eps_tau news;    // DECLARATION OF THE STRUCTURAL INNOVATIONS.

parameters phi beta kappa phi_x phi_pie alpha rho lambda; // DECLARATION OF THE PARAMETERS.

// Parameter Values  
beta    = 0.97;
phi_x   =  .0;
phi_pie =  1.5;
alpha   =  0.0;
rho     =  0.2;
lambda  =  0.5;
phi     =  1;
theta   =  0.75;

kappa   = ((1-theta)*(1-beta*theta)*(1+phi))/(theta);

// DECLARATION OF THE (LINEAR) DSGE MODEL: 
model(linear);

    beta*pie(+1) + kappa*x = pie;                                      // Calvo Pricing Equation
    r - pie(+1)-rstar   = x(+1) - x;                                   // Intertemporal Equation
    alpha*r(-1)+(1-alpha)*phi_pie*pie + (1-alpha)*phi_x*x = r;       // Monetary Policy Rule
    //r = rstar + phi_pie*pie;
    rstar = a(+1) - a + (1-lambda)/(1+phi)*tau;                           // Definition of the Natural Rate
    a = rho*a(-1) + eps_a + news(-1);
    tau       = lambda*tau(-1) + eps_tau;

end;

shocks;
var eps_a;
stderr 1;
var eps_tau;
stderr 1;
var news;
stderr 1;
end;

stoch_simul(irf=7,nograph) ;
fprintf(' r = %5.3f, rstar = %5.3f, infl = %5.3f, empl = %5.3f\n', ...
    r_news(1), rstar_news(1), pie_news(1), x_news(1));




break
tech = a_eps_a;
output = x_eps_a + tech;
tt = 1:length(x_eps_a);
ia=2;
ib=2;
subplot(ia,ib,1)
plot(tt,tech,'*-',tt,output);
legend('natural output','actual output')
axis tight;

subplot(ia,ib,2)
plot(tt,rstar_eps_a,'*-',tt,r_eps_a)
legend('natural rate', 'actual rate')
axis tight

subplot(ia,ib,3)
plot(tt,x_eps_a)
legend('employment (outout gap)')
axis tight

subplot(ia,ib,4)
plot(tt,pie_eps_a)
legend('inflation')
axis tight
