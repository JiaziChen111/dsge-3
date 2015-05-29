%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'cf';
%
% Some global variables initialization
%
global_initialization;
diary off;
diary('cf.log');
M_.exo_names = 'e';
M_.exo_names_tex = 'e';
M_.exo_names_long = 'e';
M_.exo_names = char(M_.exo_names, 's');
M_.exo_names_tex = char(M_.exo_names_tex, 's');
M_.exo_names_long = char(M_.exo_names_long, 's');
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names_long = 'y';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names_long = char(M_.endo_names_long, 'c');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names_long = char(M_.endo_names_long, 'k');
M_.endo_names = char(M_.endo_names, 'i');
M_.endo_names_tex = char(M_.endo_names_tex, 'i');
M_.endo_names_long = char(M_.endo_names_long, 'i');
M_.param_names = 'a';
M_.param_names_tex = 'a';
M_.param_names_long = 'a';
M_.param_names = char(M_.param_names, 'd');
M_.param_names_tex = char(M_.param_names_tex, 'd');
M_.param_names_long = char(M_.param_names_long, 'd');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 4;
M_.param_nbr = 2;
M_.orig_endo_nbr = 4;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.Correlation_matrix = eye(2, 2);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('cf_static');
erase_compiled_function('cf_dynamic');
M_.lead_lag_incidence = [
 0 2;
 0 3;
 1 4;
 0 5;]';
M_.nstatic = 3;
M_.nfwrd   = 0;
M_.npred   = 1;
M_.nboth   = 0;
M_.nsfwrd   = 0;
M_.nspred   = 1;
M_.ndynamic   = 1;
M_.equations_tags = {
};
M_.static_and_dynamic_models_differ = 0;
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 0;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 0;
oo_.steady_state = zeros(4, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(2, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 12;
M_.NNZDerivatives(2) = -1;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.5;
a = M_.params( 1 );
M_.params( 2 ) = 0.9;
d = M_.params( 2 );
%
% INITVAL instructions
%
options_.initval_file = 0;
oo_.steady_state( 3 ) = 0.308641975308642;
oo_.steady_state( 1 ) = 0.555555555555556;
oo_.steady_state( 4 ) = 0.277777777777778;
oo_.steady_state( 2 ) = 0.277777777777778;
oo_.exo_steady_state( 1 ) = 0;
oo_.exo_steady_state( 2 ) = 0.5;
if M_.exo_nbr > 0;
	oo_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];
end;
if M_.exo_det_nbr > 0;
	oo_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];
end;
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.01)^2;
M_.Sigma_e(2, 2) = (0.01)^2;
options_.order = 1;
options_.periods = 1000;
var_list_=[];
info = stoch_simul(var_list_);
constrained_vars_ = [];
constrained_paths_ = zeros(1, 2);
constrained_vars_ = 1;
constrained_paths_(1,1)=0.6;
constrained_paths_(1,2)=0.7;
options_cond_fcst_ = struct();
options_cond_fcst_.replic = 3000;
options_cond_fcst_.parameter_set = 'calibration';
options_cond_fcst_.controlled_varexo=[];
options_cond_fcst_.controlled_varexo = 'e';
imcforecast(constrained_paths_, constrained_vars_, options_cond_fcst_);
var_list_=[];
var_list_ = 'y';
var_list_ = char(var_list_, 'i');
plot_icforecast(var_list_, 10,options_);
save('cf_results.mat', 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save('cf_results.mat', 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save('cf_results.mat', 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save('cf_results.mat', 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save('cf_results.mat', 'estimation_info', '-append');
end


disp(['Total computing time : ' dynsec2hms(toc) ]);
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
diary off
