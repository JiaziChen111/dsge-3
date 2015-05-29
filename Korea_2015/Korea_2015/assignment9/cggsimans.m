%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'cggsimans';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('cggsimans.log');
M_.exo_names = 'eps_a';
M_.exo_names_tex = 'eps\_a';
M_.exo_names_long = 'eps_a';
M_.exo_names = char(M_.exo_names, 'eps_tau');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_tau');
M_.exo_names_long = char(M_.exo_names_long, 'eps_tau');
M_.endo_names = 'da';
M_.endo_names_tex = 'da';
M_.endo_names_long = 'da';
M_.endo_names = char(M_.endo_names, 'pie');
M_.endo_names_tex = char(M_.endo_names_tex, 'pie');
M_.endo_names_long = char(M_.endo_names_long, 'pie');
M_.endo_names = char(M_.endo_names, 'r');
M_.endo_names_tex = char(M_.endo_names_tex, 'r');
M_.endo_names_long = char(M_.endo_names_long, 'r');
M_.endo_names = char(M_.endo_names, 'rstar');
M_.endo_names_tex = char(M_.endo_names_tex, 'rstar');
M_.endo_names_long = char(M_.endo_names_long, 'rstar');
M_.endo_names = char(M_.endo_names, 'tau');
M_.endo_names_tex = char(M_.endo_names_tex, 'tau');
M_.endo_names_long = char(M_.endo_names_long, 'tau');
M_.endo_names = char(M_.endo_names, 'x');
M_.endo_names_tex = char(M_.endo_names_tex, 'x');
M_.endo_names_long = char(M_.endo_names_long, 'x');
M_.endo_names = char(M_.endo_names, 'dy');
M_.endo_names_tex = char(M_.endo_names_tex, 'dy');
M_.endo_names_long = char(M_.endo_names_long, 'dy');
M_.param_names = 'phi';
M_.param_names_tex = 'phi';
M_.param_names_long = 'phi';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names_long = char(M_.param_names_long, 'beta');
M_.param_names = char(M_.param_names, 'kappa');
M_.param_names_tex = char(M_.param_names_tex, 'kappa');
M_.param_names_long = char(M_.param_names_long, 'kappa');
M_.param_names = char(M_.param_names, 'phi_x');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_x');
M_.param_names_long = char(M_.param_names_long, 'phi_x');
M_.param_names = char(M_.param_names, 'phi_pie');
M_.param_names_tex = char(M_.param_names_tex, 'phi\_pie');
M_.param_names_long = char(M_.param_names_long, 'phi_pie');
M_.param_names = char(M_.param_names, 'alpha');
M_.param_names_tex = char(M_.param_names_tex, 'alpha');
M_.param_names_long = char(M_.param_names_long, 'alpha');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names_long = char(M_.param_names_long, 'rho');
M_.param_names = char(M_.param_names, 'lambda');
M_.param_names_tex = char(M_.param_names_tex, 'lambda');
M_.param_names_long = char(M_.param_names_long, 'lambda');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 7;
M_.param_nbr = 8;
M_.orig_endo_nbr = 7;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('cggsimans_static');
erase_compiled_function('cggsimans_dynamic');
M_.lead_lag_incidence = [
 1 5 12;
 0 6 13;
 2 7 0;
 0 8 0;
 3 9 0;
 4 10 14;
 0 11 0;]';
M_.nstatic = 2;
M_.nfwrd   = 1;
M_.npred   = 2;
M_.nboth   = 2;
M_.nsfwrd   = 3;
M_.nspred   = 4;
M_.ndynamic   = 5;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(7, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(8, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 27;
M_.NNZDerivatives(2) = 0;
M_.NNZDerivatives(3) = -1;
M_.params( 2 ) = 0.97;
beta = M_.params( 2 );
M_.params( 4 ) = .15;
phi_x = M_.params( 4 );
M_.params( 5 ) = 1.5;
phi_pie = M_.params( 5 );
M_.params( 6 ) = 0.8;
alpha = M_.params( 6 );
M_.params( 7 ) = 0.9;
rho = M_.params( 7 );
M_.params( 8 ) = 0.5;
lambda = M_.params( 8 );
M_.params( 1 ) = 1;
phi = M_.params( 1 );
theta   = 0.75;
M_.params( 3 ) = (1-theta)*(1-theta*M_.params(2))*(1+M_.params(1))/theta;
kappa = M_.params( 3 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (.02)^2;
M_.Sigma_e(2, 2) = (.02)^2;
set_dynare_seed=1;
options_.irf = 7;
options_.nograph = 1;
options_.periods = 5000;
var_list_=[];
info = stoch_simul(var_list_);
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
subplot(211)
lam=1600;
[y_hp,y_hptrend]=HPFAST(ysim,lam);
ystar=ysim-gap;
ca=corrcoef(gap,y_hp);
sgap=std(gap);
syhp=std(y_hp);
plot(tt,y_hptrend(tt),tt,ystar(tt),'+-',tt,ysim(tt));
legend('hp trend','natural output','actual output')
axis tight;
subplot(212)
plot(tt,y_hp(tt),tt,gap(tt),'+-')
legend('hp-filtered output','actual gap')
str=['correlation, hp-filtered output and actual output gap = ',num2str(ca(1,2)),' std(gap) = ',num2str(sgap),' std(yhp) = ',num2str(syhp)];
title(str)
axis tight;
sim = oo_.endo_simul;
save data sim;
save('cggsimans_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('cggsimans_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('cggsimans_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('cggsimans_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('cggsimans_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
