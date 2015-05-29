// Maximum Likelihood estimation of the Clarida, Gali, Gertler Model

var da pie r rstar tau x dy; 	   // DECLARATION OF THE ENDOGENOUS VARIABLES. 

varexo eps_a eps_tau;    // DECLARATION OF THE STRUCTURAL INNOVATIONS.

parameters phi beta kappa phi_x phi_pie alpha rho lambda; 

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
    rstar = rho*da + (1-lambda)/(1+phi)*tau;                           // Definition of the Natural Rate  
    da = rho*da(-1) + eps_a;                                           // Exogenous Stochastic Processes
    tau       = lambda*tau(-1) + eps_tau;
    dy  = x - x(-1) + da - (tau - tau(-1))/(1+phi);                    // construction of output growth

end;


// maximum likelihood estimation:
estimated_params;
   stderr eps_a, 0.02;
   stderr eps_tau, 0.02;
   rho, .9;
   lambda, .5;  
end;


estimated_params_bounds;
stderr eps_a, 0.001, .2;
stderr eps_tau, 0.001, .2;
rho, .001,.95;
lambda, .001,.95;
end;


// Names of variables observed in the estimation
varobs pie dy; 

// arguments in the call to the Dynare command, estimation.
// datafile = xx, xx is an *.m file which, when executed, produces variables, pie and dy (at least)
// conf_sig = xx, xx defines width of confidence intervals (e.g., xx = .95 means '95 percent confidence interval')
// first_obs=xx, xx is the first observation in the data vectors produced by datapullcgg, which is used in estimation (xx must be an actual number, not a pre-defined variable)
// nobs=yy, Dynare uses observations t=xx to t=yy-1 in estimation (yy must be an actual number)
// mode_check, activates a graph which plots the posterior mode against the criterion. This graph provides visual confirmation that the function of interest in fact has been optimized
// mode_compute=xx, xx selects one of several possible optimization routines to use in optimization
// forecast triggers computation of forcasts n periods after the end of the data set, were forecast=n.


estimation(datafile=datapullcgg,conf_sig =.95,first_obs=51,forecast =36,nobs=30,mode_check,mode_compute=4) da pie r rstar tau x dy;

