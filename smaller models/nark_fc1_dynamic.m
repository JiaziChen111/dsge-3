function [residual, g1, g2, g3] = nark_fc1_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(19, 1);
T54 = exp(y(12))-params(5)*exp(y(1))/exp(y(22));
T61 = exp(y(31))-exp(y(12))*params(5)/exp(y(32));
T71 = params(29)/(params(29)-1);
T76 = exp(y(15))^params(25);
T107 = exp(y(2))*exp(y(20))*exp(y(23))*(1-params(4))/params(4);
T109 = T107/exp(y(15))/exp(y(22));
T111 = exp(y(14))/(1-params(4));
T114 = (exp(y(20))/params(4))^params(4);
T121 = exp(y(2))^params(4)*exp(y(15))^(1-params(4));
T123 = exp(y(22))^(-params(4));
T164 = exp((-params(14))*(exp(y(27))*exp(y(29))-exp((steady_state(7)))*params(13))-params(33)*(exp(y(25))-exp((steady_state(14)))-(exp(y(13))-exp((steady_state(2))))));
T296 = (-((-(params(5)*exp(y(1))*exp(y(22))))/(exp(y(22))*exp(y(22)))));
lhs =y(22);
rhs =params(36)*y(4)+(1-params(36))*params(18)+x(it_, 3);
residual(1)= lhs-rhs;
lhs =y(21);
rhs =(1-params(49))*params(6)+params(49)*y(3)+x(it_, 17);
residual(2)= lhs-rhs;
lhs =y(23);
rhs =params(47)*y(5)+x(it_, 16);
residual(3)= lhs-rhs;
lhs =y(30);
rhs =params(48)*y(11)+x(it_, 15);
residual(4)= lhs-rhs;
lhs =exp(y(32))/T54;
rhs =params(1)/T61*exp(y(13));
residual(5)= lhs-rhs;
lhs =exp(y(14));
rhs =T54*T71*T76;
residual(6)= lhs-rhs;
lhs =1/exp(y(16));
rhs =y(21)/(y(21)-1);
residual(7)= lhs-rhs;
lhs =exp(y(20));
rhs =exp(y(13))-(1-params(2));
residual(8)= lhs-rhs;
lhs =exp(y(17))-(1-params(2))*exp(y(2))/exp(y(22));
rhs =exp(y(19));
residual(9)= lhs-rhs;
lhs =exp(y(14));
rhs =T109;
residual(10)= lhs-rhs;
lhs =exp(y(16));
rhs =T111^(1-params(4))*T114;
residual(11)= lhs-rhs;
lhs =exp(y(18));
rhs =T121*T123;
residual(12)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(12))+exp(y(19));
residual(13)= lhs-rhs;
lhs =(exp(y(20))+1-params(2))*exp(y(26))*exp(y(27));
rhs =exp(y(28))*exp(y(33));
residual(14)= lhs-rhs;
lhs =exp(y(28));
rhs =exp(y(25))*T164*exp(y(30));
residual(15)= lhs-rhs;
lhs =exp(y(27))*exp(y(9))*exp(y(10));
rhs =exp(y(26))*exp(y(22))*exp(y(27))*exp(y(29));
residual(16)= lhs-rhs;
lhs =y(24);
rhs =params(61)*y(6)+params(62)*y(7)+params(63)*y(8)+x(it_, 12);
residual(17)= lhs-rhs;
lhs =y(25);
rhs =y(6)*params(64)+y(7)*params(65)+y(8)*params(66)+x(it_, 13);
residual(18)= lhs-rhs;
lhs =y(26);
rhs =y(6)*params(67)+y(7)*params(68)+y(8)*params(69)+x(it_, 14);
residual(19)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(19, 51);

  %
  % Jacobian matrix
  %

  g1(1,4)=(-params(36));
  g1(1,22)=1;
  g1(1,36)=(-1);
  g1(2,3)=(-params(49));
  g1(2,21)=1;
  g1(2,50)=(-1);
  g1(3,5)=(-params(47));
  g1(3,23)=1;
  g1(3,49)=(-1);
  g1(4,11)=(-params(48));
  g1(4,30)=1;
  g1(4,48)=(-1);
  g1(5,1)=(-(exp(y(32))*(-(params(5)*exp(y(1))/exp(y(22))))))/(T54*T54);
  g1(5,12)=(-(exp(y(32))*exp(y(12))))/(T54*T54)-exp(y(13))*(-(params(1)*(-(exp(y(12))*params(5)/exp(y(32))))))/(T61*T61);
  g1(5,31)=(-(exp(y(13))*(-(params(1)*exp(y(31))))/(T61*T61)));
  g1(5,13)=(-(params(1)/T61*exp(y(13))));
  g1(5,22)=(-(exp(y(32))*T296))/(T54*T54);
  g1(5,32)=exp(y(32))/T54-exp(y(13))*(-(params(1)*(-((-(exp(y(32))*exp(y(12))*params(5)))/(exp(y(32))*exp(y(32)))))))/(T61*T61);
  g1(6,1)=(-(T76*T71*(-(params(5)*exp(y(1))/exp(y(22))))));
  g1(6,12)=(-(T76*exp(y(12))*T71));
  g1(6,14)=exp(y(14));
  g1(6,15)=(-(T54*T71*exp(y(15))*getPowerDeriv(exp(y(15)),params(25),1)));
  g1(6,22)=(-(T76*T71*T296));
  g1(7,16)=(-exp(y(16)))/(exp(y(16))*exp(y(16)));
  g1(7,21)=(-((y(21)-1-y(21))/((y(21)-1)*(y(21)-1))));
  g1(8,13)=(-exp(y(13)));
  g1(8,20)=exp(y(20));
  g1(9,2)=(-((1-params(2))*exp(y(2))/exp(y(22))));
  g1(9,17)=exp(y(17));
  g1(9,19)=(-exp(y(19)));
  g1(9,22)=(-((-(exp(y(22))*(1-params(2))*exp(y(2))))/(exp(y(22))*exp(y(22)))));
  g1(10,14)=exp(y(14));
  g1(10,15)=(-((-(exp(y(15))*T107))/(exp(y(15))*exp(y(15)))/exp(y(22))));
  g1(10,2)=(-T109);
  g1(10,20)=(-T109);
  g1(10,22)=(-((-(exp(y(22))*T107/exp(y(15))))/(exp(y(22))*exp(y(22)))));
  g1(10,23)=(-T109);
  g1(11,14)=(-(T114*T111*getPowerDeriv(T111,1-params(4),1)));
  g1(11,16)=exp(y(16));
  g1(11,20)=(-(T111^(1-params(4))*exp(y(20))/params(4)*getPowerDeriv(exp(y(20))/params(4),params(4),1)));
  g1(12,15)=(-(T123*exp(y(2))^params(4)*exp(y(15))*getPowerDeriv(exp(y(15)),1-params(4),1)));
  g1(12,2)=(-(T123*exp(y(15))^(1-params(4))*exp(y(2))*getPowerDeriv(exp(y(2)),params(4),1)));
  g1(12,18)=exp(y(18));
  g1(12,22)=(-(T121*exp(y(22))*getPowerDeriv(exp(y(22)),(-params(4)),1)));
  g1(13,12)=(-exp(y(12)));
  g1(13,18)=exp(y(18));
  g1(13,19)=(-exp(y(19)));
  g1(14,20)=exp(y(27))*exp(y(20))*exp(y(26));
  g1(14,26)=(exp(y(20))+1-params(2))*exp(y(26))*exp(y(27));
  g1(14,27)=(exp(y(20))+1-params(2))*exp(y(26))*exp(y(27));
  g1(14,33)=(-(exp(y(28))*exp(y(33))));
  g1(14,28)=(-(exp(y(28))*exp(y(33))));
  g1(15,13)=(-(exp(y(25))*exp(y(30))*T164*(-(params(33)*(-exp(y(13)))))));
  g1(15,25)=(-(exp(y(25))*T164*exp(y(30))+exp(y(25))*exp(y(30))*T164*(-(params(33)*exp(y(25))))));
  g1(15,27)=(-(exp(y(25))*exp(y(30))*T164*(-params(14))*exp(y(27))*exp(y(29))));
  g1(15,28)=exp(y(28));
  g1(15,29)=(-(exp(y(25))*exp(y(30))*T164*(-params(14))*exp(y(27))*exp(y(29))));
  g1(15,30)=(-(exp(y(25))*T164*exp(y(30))));
  g1(16,22)=(-(exp(y(26))*exp(y(22))*exp(y(27))*exp(y(29))));
  g1(16,26)=(-(exp(y(26))*exp(y(22))*exp(y(27))*exp(y(29))));
  g1(16,27)=exp(y(27))*exp(y(9))*exp(y(10))-exp(y(26))*exp(y(22))*exp(y(27))*exp(y(29));
  g1(16,9)=exp(y(27))*exp(y(9))*exp(y(10));
  g1(16,10)=exp(y(27))*exp(y(9))*exp(y(10));
  g1(16,29)=(-(exp(y(26))*exp(y(22))*exp(y(27))*exp(y(29))));
  g1(17,6)=(-params(61));
  g1(17,24)=1;
  g1(17,7)=(-params(62));
  g1(17,8)=(-params(63));
  g1(17,45)=(-1);
  g1(18,6)=(-params(64));
  g1(18,7)=(-params(65));
  g1(18,25)=1;
  g1(18,8)=(-params(66));
  g1(18,46)=(-1);
  g1(19,6)=(-params(67));
  g1(19,7)=(-params(68));
  g1(19,8)=(-params(69));
  g1(19,26)=1;
  g1(19,47)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],19,2601);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],19,132651);
end
end
