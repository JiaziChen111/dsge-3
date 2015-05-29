// Estimation of the Clarida, Gali, Gertler Model

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


// Set up the priors on the parameters
// Bayesian estimation:
estimated_params;
   stderr eps_a, inv_gamma_pdf,0.1,1;
   stderr eps_tau, inv_gamma_pdf,0.1,1;
   rho, beta_pdf, 0.9, .04;
   lambda, beta_pdf, 0.5, .04; 
end;

// maximum likelihood estimation:
//estimated_params;
//   stderr eps_a, 0.02;
//   stderr eps_tau, 0.02;
//   rho, .95;
//   lambda, .5;  
//end;

estimated_params_bounds;
stderr eps_a, 0.001, .2;
stderr eps_tau, 0.001, .2;
rho, .001,.99;
lambda, .001,.99;
end;

// the following is optional
// Initial values for the parameters in the estimation, pre-set to the optimum..
estimated_params_init;
stderr eps_a, 0.02;
stderr eps_tau, 0.02;
rho, 0.9;
lambda, 0.5;
end;

// Names of variables observed in the estimation
varobs pie dy; 

// arguments in the call to the Dynare command, estimation.
// datafile = xx, xx is an *.m file which, when executed, produces variables, pie and dy (at least)
// conf_sig = xx, xx defines width of confidence intervals
// smoother, triggers computation of 'forecast' of state at each date, based on all observed variables, past and future. 
// first_obs=xx, xx is the first observation in the data vectors produced by datapullcgg, which is used in estimation (xx must be an actual number, not a pre-defined variable)
// nobs=yy, Dynare uses observations t=xx to t=xx+yy-1 in estimation (yy must be an actual number)
// mode_check, activates a graph which plots the posterior mode against the criterion. This graph provides visual confirmation that the function of interest in fact has been optimized
// mode_compute=xx, xx selects one of several possible optimization routines to use in optimization
// mh_replic=xx, xx is the number of Metropolis-Hastings trials to do in computing the posterior distribution
// mh_jscale=xx, xx is the scale factor to use in Metropolis-Hastings trials. The 'acceptation rate' should be around 30%. If it is higher, raise mh_jscale
// mh_nblocks=xx, xx is the number of times to repeat the Metropolis-Hastings trials.
// filtered_vars, triggers computation of 'forecast' of state at each date, based on all previous observations only.
// Triggers the computation of the posterior distribution of the theoretical moments of the endogenous variables. Results are stored in oo_.PosteriorTheoreticalMoments
// filtered_vars triggers the computation of the posterior distribution of filtered endogenous variables/one-step ahead forecasts. Results are stored in oo_.FilteredVariables.
// forecast triggers computation of forcasts n periods after the end of the data set, were forecast=n.


// in case you run analyzegap, make sure the values here correspond to first_obs and nobs in the estimation command:
first_=51;
%num_=3000;
// Bayesian estimation 
estimation(datafile=datapullcgg,conf_sig =.95,smoother,moments_varendo,first_obs=51,forecast =36,nobs=30,mode_check,mode_compute=4,mh_replic=1000,mh_jscale=1.5,mh_nblocks=1) da pie r rstar tau x dy;

//the following lines should only be run if you've done Bayesian estimation
compareMCMCLaplace
//figure
// evaluate alternative estimates of the output gap
//analyzegap;
