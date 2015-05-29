function [residual, g1, g2, g3] = nark_fc_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [M_.exo_nbr by nperiods] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(14, 1);
T47 = exp(y(9))-params(5)*exp(y(1))/exp(y(18));
T54 = exp(y(23))-exp(y(9))*params(5)/exp(y(24));
T67 = params(29)/(params(29)-1);
T72 = exp(y(11))^params(25);
T97 = exp(y(2))*exp(y(16))*exp(y(19))*(1-params(4))/params(4);
T99 = T97/exp(y(11))/exp(y(18));
T101 = exp(y(10))/(1-params(4));
T104 = (exp(y(16))/params(4))^params(4);
T111 = exp(y(2))^params(4)*exp(y(11))^(1-params(4));
T113 = exp(y(18))^(-params(4));
T230 = (-((-(params(5)*exp(y(1))*exp(y(18))))/(exp(y(18))*exp(y(18)))));
lhs =y(18);
rhs =params(36)*y(4)+(1-params(36))*params(18)+x(it_, 3);
residual(1)= lhs-rhs;
lhs =y(17);
rhs =(1-params(49))*params(6)+params(49)*y(3)+x(it_, 17);
residual(2)= lhs-rhs;
lhs =y(19);
rhs =params(47)*y(5)+x(it_, 16);
residual(3)= lhs-rhs;
lhs =exp(y(24))/T47;
rhs =params(1)/T54*(exp(y(16))+1-params(2));
residual(4)= lhs-rhs;
lhs =exp(y(10));
rhs =T47*T67*T72;
residual(5)= lhs-rhs;
lhs =1/exp(y(12));
rhs =y(17)/(y(17)-1);
residual(6)= lhs-rhs;
lhs =exp(y(13))-(1-params(2))*exp(y(2))/exp(y(18));
rhs =exp(y(15));
residual(7)= lhs-rhs;
lhs =exp(y(10));
rhs =T99;
residual(8)= lhs-rhs;
lhs =exp(y(12));
rhs =T101^(1-params(4))*T104;
residual(9)= lhs-rhs;
lhs =exp(y(14));
rhs =T111*T113;
residual(10)= lhs-rhs;
lhs =exp(y(14));
rhs =exp(y(9))+exp(y(15));
residual(11)= lhs-rhs;
lhs =y(20);
rhs =params(61)*y(6)+params(62)*y(7)+params(63)*y(8)+x(it_, 12);
residual(12)= lhs-rhs;
lhs =y(21);
rhs =y(6)*params(64)+y(7)*params(65)+y(8)*params(66)+x(it_, 13);
residual(13)= lhs-rhs;
lhs =y(22);
rhs =y(6)*params(67)+y(7)*params(68)+y(8)*params(69)+x(it_, 14);
residual(14)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(14, 42);

  %
  % Jacobian matrix
  %

  g1(1,4)=(-params(36));
  g1(1,18)=1;
  g1(1,27)=(-1);
  g1(2,3)=(-params(49));
  g1(2,17)=1;
  g1(2,41)=(-1);
  g1(3,5)=(-params(47));
  g1(3,19)=1;
  g1(3,40)=(-1);
  g1(4,1)=(-(exp(y(24))*(-(params(5)*exp(y(1))/exp(y(18))))))/(T47*T47);
  g1(4,9)=(-(exp(y(24))*exp(y(9))))/(T47*T47)-(exp(y(16))+1-params(2))*(-(params(1)*(-(exp(y(9))*params(5)/exp(y(24))))))/(T54*T54);
  g1(4,23)=(-((exp(y(16))+1-params(2))*(-(params(1)*exp(y(23))))/(T54*T54)));
  g1(4,16)=(-(params(1)/T54*exp(y(16))));
  g1(4,18)=(-(exp(y(24))*T230))/(T47*T47);
  g1(4,24)=exp(y(24))/T47-(exp(y(16))+1-params(2))*(-(params(1)*(-((-(exp(y(24))*exp(y(9))*params(5)))/(exp(y(24))*exp(y(24)))))))/(T54*T54);
  g1(5,1)=(-(T72*T67*(-(params(5)*exp(y(1))/exp(y(18))))));
  g1(5,9)=(-(T72*exp(y(9))*T67));
  g1(5,10)=exp(y(10));
  g1(5,11)=(-(T47*T67*exp(y(11))*getPowerDeriv(exp(y(11)),params(25),1)));
  g1(5,18)=(-(T72*T67*T230));
  g1(6,12)=(-exp(y(12)))/(exp(y(12))*exp(y(12)));
  g1(6,17)=(-((y(17)-1-y(17))/((y(17)-1)*(y(17)-1))));
  g1(7,2)=(-((1-params(2))*exp(y(2))/exp(y(18))));
  g1(7,13)=exp(y(13));
  g1(7,15)=(-exp(y(15)));
  g1(7,18)=(-((-(exp(y(18))*(1-params(2))*exp(y(2))))/(exp(y(18))*exp(y(18)))));
  g1(8,10)=exp(y(10));
  g1(8,11)=(-((-(exp(y(11))*T97))/(exp(y(11))*exp(y(11)))/exp(y(18))));
  g1(8,2)=(-T99);
  g1(8,16)=(-T99);
  g1(8,18)=(-((-(exp(y(18))*T97/exp(y(11))))/(exp(y(18))*exp(y(18)))));
  g1(8,19)=(-T99);
  g1(9,10)=(-(T104*T101*getPowerDeriv(T101,1-params(4),1)));
  g1(9,12)=exp(y(12));
  g1(9,16)=(-(T101^(1-params(4))*exp(y(16))/params(4)*getPowerDeriv(exp(y(16))/params(4),params(4),1)));
  g1(10,11)=(-(T113*exp(y(2))^params(4)*exp(y(11))*getPowerDeriv(exp(y(11)),1-params(4),1)));
  g1(10,2)=(-(T113*exp(y(11))^(1-params(4))*exp(y(2))*getPowerDeriv(exp(y(2)),params(4),1)));
  g1(10,14)=exp(y(14));
  g1(10,18)=(-(T111*exp(y(18))*getPowerDeriv(exp(y(18)),(-params(4)),1)));
  g1(11,9)=(-exp(y(9)));
  g1(11,14)=exp(y(14));
  g1(11,15)=(-exp(y(15)));
  g1(12,6)=(-params(61));
  g1(12,20)=1;
  g1(12,7)=(-params(62));
  g1(12,8)=(-params(63));
  g1(12,36)=(-1);
  g1(13,6)=(-params(64));
  g1(13,7)=(-params(65));
  g1(13,21)=1;
  g1(13,8)=(-params(66));
  g1(13,37)=(-1);
  g1(14,6)=(-params(67));
  g1(14,7)=(-params(68));
  g1(14,8)=(-params(69));
  g1(14,22)=1;
  g1(14,38)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],14,1764);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],14,74088);
end
end
