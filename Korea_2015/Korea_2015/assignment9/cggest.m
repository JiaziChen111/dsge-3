%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'cggest';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('cggest.log');
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
options_.varobs = [];
options_.varobs = 'pie';
options_.varobs = char(options_.varobs, 'dy');
options_.varobs_id = [ 2 7  ];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('cggest_static');
erase_compiled_function('cggest_dynamic');
M_.lead_lag_incidence = [
 1 5 0;
 0 6 12;
 2 7 0;
 0 8 0;
 3 9 0;
 4 10 13;
 0 11 0;]';
M_.nstatic = 2;
M_.nfwrd   = 1;
M_.npred   = 3;
M_.nboth   = 1;
M_.nsfwrd   = 2;
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
global estim_params_
estim_params_.var_exo = [];
estim_params_.var_endo = [];
estim_params_.corrx = [];
estim_params_.corrn = [];
estim_params_.param_vals = [];
estim_params_.var_exo = [estim_params_.var_exo; 1, 0.02, (-Inf), Inf, 0, NaN, NaN, NaN, NaN, NaN ];
estim_params_.var_exo = [estim_params_.var_exo; 2, 0.02, (-Inf), Inf, 0, NaN, NaN, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 7, .9, (-Inf), Inf, 0, NaN, NaN, NaN, NaN, NaN ];
estim_params_.param_vals = [estim_params_.param_vals; 8, .5, (-Inf), Inf, 0, NaN, NaN, NaN, NaN, NaN ];
tmp1 = find(estim_params_.var_exo(:,1)==1);
estim_params_.var_exo(tmp1,3) = 0.001;
estim_params_.var_exo(tmp1,4) = .2;
tmp1 = find(estim_params_.var_exo(:,1)==2);
estim_params_.var_exo(tmp1,3) = 0.001;
estim_params_.var_exo(tmp1,4) = .2;
tmp1 = find(estim_params_.param_vals(:,1)==7);
estim_params_.param_vals(tmp1,3) = .001;
estim_params_.param_vals(tmp1,4) = .95;
tmp1 = find(estim_params_.param_vals(:,1)==8);
estim_params_.param_vals(tmp1,3) = .001;
estim_params_.param_vals(tmp1,4) = .95;
options_.conf_sig = .95;
options_.first_obs = 51;
options_.forecast = 36;
options_.mode_check.status = 1;
options_.mode_compute = 4;
options_.datafile = 'datapullcgg';
options_.nobs = 30;
options_.order = 1;
var_list_=[];
var_list_ = 'da';
var_list_ = char(var_list_, 'pie');
var_list_ = char(var_list_, 'r');
var_list_ = char(var_list_, 'rstar');
var_list_ = char(var_list_, 'tau');
var_list_ = char(var_list_, 'x');
var_list_ = char(var_list_, 'dy');
dynare_estimation(var_list_);
save('cggest_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('cggest_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('cggest_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('cggest_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('cggest_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
