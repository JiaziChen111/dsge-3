function [residual, g1, g2, g3] = nark_fc_fb_xm_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(34, 1);
T19 = exp(y(14))-params(5)*exp(y(1))/exp(y(24));
T26 = exp(y(48))-exp(y(14))*params(5)/exp(y(49));
T36 = params(29)/(params(29)-1);
T41 = exp(y(17))^params(25);
T67 = exp(y(2))*exp(y(22))*exp(y(25))*(1-params(4))/params(4);
T69 = T67/exp(y(17))/exp(y(24));
T73 = exp(y(16))/(1-params(4));
T76 = (exp(y(22))/params(4))^params(4);
T83 = exp(y(2))^params(4)*exp(y(17))^(1-params(4));
T85 = exp(y(24))^(-params(4));
T133 = exp((-params(14))*(exp(y(29))*exp(y(31))-exp((steady_state(7)))*params(13))-params(33)*(exp(y(27))-exp((steady_state(14)))-(exp(y(15))-exp((steady_state(2))))));
T150 = exp(y(29))*exp(y(9))*exp(y(10))/exp(y(24));
T151 = T150/exp(y(28));
T180 = exp(y(14))*(1-params(21))*exp(y(36))^(-params(22));
T197 = exp(y(37))*params(7)/(params(7)-1);
T204 = exp(y(12))/exp(y(24));
T207 = exp(y(47))*T204^params(15);
T213 = exp(y(45))*(1+params(30)*(exp(y(15))-1))/exp(y(29));
T219 = T213^(-params(16))*exp(y(26));
T221 = T219^(1-params(15));
T225 = exp(y(35))/exp(y(45));
T235 = exp(y(46))/exp(y(45));
T238 = exp(y(42))*(1-params(34))*T235^(-params(35));
T245 = params(34)*exp(y(35))^(1-params(35))+(1-params(34))*exp(y(46))^(1-params(35));
T253 = exp(y(35))/exp(y(40));
T258 = exp(y(21))*params(23)*T253^(-params(24));
T263 = exp(y(41))/exp(y(40));
T266 = exp(y(21))*(1-params(23))*T263^(-params(24));
T273 = params(23)*exp(y(35))^(1-params(24))+(1-params(23))*exp(y(41))^(1-params(24));
T391 = getPowerDeriv(T213,(-params(16)),1);
T394 = getPowerDeriv(T219,1-params(15),1);
T445 = (-((-(params(5)*exp(y(1))*exp(y(24))))/(exp(y(24))*exp(y(24)))));
T550 = getPowerDeriv(T245,1/(1-params(35)),1);
T561 = getPowerDeriv(T273,1/(1-params(24)),1);
lhs =exp(y(49))/T19;
rhs =params(1)/T26*exp(y(15));
residual(1)= lhs-rhs;
lhs =exp(y(16));
rhs =T19*T36*T41;
residual(2)= lhs-rhs;
lhs =exp(y(22));
rhs =exp(y(15))-(1-params(2));
residual(3)= lhs-rhs;
lhs =exp(y(19))-(1-params(2))*exp(y(2))/exp(y(24));
rhs =exp(y(21));
residual(4)= lhs-rhs;
lhs =exp(y(16));
rhs =T69;
residual(5)= lhs-rhs;
lhs =exp(y(18));
rhs =T73^(1-params(4))*T76;
residual(6)= lhs-rhs;
lhs =exp(y(20));
rhs =T83*T85;
residual(7)= lhs-rhs;
lhs =exp(y(20));
rhs =exp(y(33))+exp(y(38))+exp(y(43));
residual(8)= lhs-rhs;
lhs =(exp(y(22))+1-params(2))*exp(y(28))*exp(y(29));
rhs =exp(y(30))*exp(y(50));
residual(9)= lhs-rhs;
lhs =exp(y(30));
rhs =exp(y(27))*T133*exp(y(32));
residual(10)= lhs-rhs;
lhs =exp(y(45))*exp(y(42))+T151;
rhs =exp(y(29))*exp(y(31))+exp(y(37))*(exp(y(34))+exp(y(44))+exp(y(39)));
residual(11)= lhs-rhs;
lhs =exp(y(33));
rhs =exp(y(14))*params(21)*exp(y(35))^(-params(22));
residual(12)= lhs-rhs;
lhs =exp(y(34));
rhs =T180;
residual(13)= lhs-rhs;
lhs =1;
rhs =params(21)*exp(y(35))^(1-params(22))+(1-params(21))*exp(y(36))^(1-params(22));
residual(14)= lhs-rhs;
lhs =exp(y(35));
rhs =exp(y(18))*y(23)/(y(23)-1);
residual(15)= lhs-rhs;
lhs =exp(y(36));
rhs =T197;
residual(16)= lhs-rhs;
lhs =exp(y(37));
rhs =exp(y(29));
residual(17)= lhs-rhs;
lhs =exp(y(42));
rhs =T207*T221;
residual(18)= lhs-rhs;
lhs =exp(y(43));
rhs =exp(y(42))*params(34)*T225^(-params(35));
residual(19)= lhs-rhs;
lhs =exp(y(44));
rhs =T238;
residual(20)= lhs-rhs;
lhs =exp(y(45));
rhs =T245^(1/(1-params(35)));
residual(21)= lhs-rhs;
lhs =exp(y(46));
rhs =exp(y(37));
residual(22)= lhs-rhs;
lhs =exp(y(38));
rhs =T258;
residual(23)= lhs-rhs;
lhs =exp(y(39));
rhs =T266;
residual(24)= lhs-rhs;
lhs =exp(y(40));
rhs =T273^(1/(1-params(24)));
residual(25)= lhs-rhs;
lhs =exp(y(41));
rhs =T197;
residual(26)= lhs-rhs;
lhs =y(26);
rhs =params(61)*y(6)+params(62)*y(7)+params(63)*y(8)+x(it_, 12);
residual(27)= lhs-rhs;
lhs =y(27);
rhs =y(6)*params(64)+y(7)*params(65)+y(8)*params(66)+x(it_, 13);
residual(28)= lhs-rhs;
lhs =y(28);
rhs =y(6)*params(67)+y(7)*params(68)+y(8)*params(69)+x(it_, 14);
residual(29)= lhs-rhs;
lhs =y(24);
rhs =params(36)*y(4)+(1-params(36))*params(18)+x(it_, 3);
residual(30)= lhs-rhs;
lhs =y(23);
rhs =(1-params(49))*params(6)+params(49)*y(3)+x(it_, 17);
residual(31)= lhs-rhs;
lhs =y(25);
rhs =params(47)*y(5)+x(it_, 16);
residual(32)= lhs-rhs;
lhs =y(32);
rhs =params(48)*y(11)+x(it_, 15);
residual(33)= lhs-rhs;
lhs =y(47);
rhs =params(43)*y(13)+x(it_, 9);
residual(34)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(34, 68);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-(exp(y(49))*(-(params(5)*exp(y(1))/exp(y(24))))))/(T19*T19);
  g1(1,14)=(-(exp(y(49))*exp(y(14))))/(T19*T19)-exp(y(15))*(-(params(1)*(-(exp(y(14))*params(5)/exp(y(49))))))/(T26*T26);
  g1(1,48)=(-(exp(y(15))*(-(params(1)*exp(y(48))))/(T26*T26)));
  g1(1,15)=(-(params(1)/T26*exp(y(15))));
  g1(1,24)=(-(exp(y(49))*T445))/(T19*T19);
  g1(1,49)=exp(y(49))/T19-exp(y(15))*(-(params(1)*(-((-(exp(y(49))*exp(y(14))*params(5)))/(exp(y(49))*exp(y(49)))))))/(T26*T26);
  g1(2,1)=(-(T41*T36*(-(params(5)*exp(y(1))/exp(y(24))))));
  g1(2,14)=(-(T41*exp(y(14))*T36));
  g1(2,16)=exp(y(16));
  g1(2,17)=(-(T19*T36*exp(y(17))*getPowerDeriv(exp(y(17)),params(25),1)));
  g1(2,24)=(-(T41*T36*T445));
  g1(3,15)=(-exp(y(15)));
  g1(3,22)=exp(y(22));
  g1(4,2)=(-((1-params(2))*exp(y(2))/exp(y(24))));
  g1(4,19)=exp(y(19));
  g1(4,21)=(-exp(y(21)));
  g1(4,24)=(-((-(exp(y(24))*(1-params(2))*exp(y(2))))/(exp(y(24))*exp(y(24)))));
  g1(5,16)=exp(y(16));
  g1(5,17)=(-((-(exp(y(17))*T67))/(exp(y(17))*exp(y(17)))/exp(y(24))));
  g1(5,2)=(-T69);
  g1(5,22)=(-T69);
  g1(5,24)=(-((-(exp(y(24))*T67/exp(y(17))))/(exp(y(24))*exp(y(24)))));
  g1(5,25)=(-T69);
  g1(6,16)=(-(T76*T73*getPowerDeriv(T73,1-params(4),1)));
  g1(6,18)=exp(y(18));
  g1(6,22)=(-(T73^(1-params(4))*exp(y(22))/params(4)*getPowerDeriv(exp(y(22))/params(4),params(4),1)));
  g1(7,17)=(-(T85*exp(y(2))^params(4)*exp(y(17))*getPowerDeriv(exp(y(17)),1-params(4),1)));
  g1(7,2)=(-(T85*exp(y(17))^(1-params(4))*exp(y(2))*getPowerDeriv(exp(y(2)),params(4),1)));
  g1(7,20)=exp(y(20));
  g1(7,24)=(-(T83*exp(y(24))*getPowerDeriv(exp(y(24)),(-params(4)),1)));
  g1(8,20)=exp(y(20));
  g1(8,33)=(-exp(y(33)));
  g1(8,38)=(-exp(y(38)));
  g1(8,43)=(-exp(y(43)));
  g1(9,22)=exp(y(29))*exp(y(22))*exp(y(28));
  g1(9,28)=(exp(y(22))+1-params(2))*exp(y(28))*exp(y(29));
  g1(9,29)=(exp(y(22))+1-params(2))*exp(y(28))*exp(y(29));
  g1(9,50)=(-(exp(y(30))*exp(y(50))));
  g1(9,30)=(-(exp(y(30))*exp(y(50))));
  g1(10,15)=(-(exp(y(27))*exp(y(32))*T133*(-(params(33)*(-exp(y(15)))))));
  g1(10,27)=(-(exp(y(27))*T133*exp(y(32))+exp(y(27))*exp(y(32))*T133*(-(params(33)*exp(y(27))))));
  g1(10,29)=(-(exp(y(27))*exp(y(32))*T133*(-params(14))*exp(y(29))*exp(y(31))));
  g1(10,30)=exp(y(30));
  g1(10,31)=(-(exp(y(27))*exp(y(32))*T133*(-params(14))*exp(y(29))*exp(y(31))));
  g1(10,32)=(-(exp(y(27))*T133*exp(y(32))));
  g1(11,24)=(-(exp(y(24))*exp(y(29))*exp(y(9))*exp(y(10))))/(exp(y(24))*exp(y(24)))/exp(y(28));
  g1(11,28)=(-(exp(y(28))*T150))/(exp(y(28))*exp(y(28)));
  g1(11,29)=T151-exp(y(29))*exp(y(31));
  g1(11,9)=T151;
  g1(11,10)=T151;
  g1(11,31)=(-(exp(y(29))*exp(y(31))));
  g1(11,34)=(-(exp(y(37))*exp(y(34))));
  g1(11,37)=(-(exp(y(37))*(exp(y(34))+exp(y(44))+exp(y(39)))));
  g1(11,39)=(-(exp(y(37))*exp(y(39))));
  g1(11,42)=exp(y(45))*exp(y(42));
  g1(11,44)=(-(exp(y(37))*exp(y(44))));
  g1(11,45)=exp(y(45))*exp(y(42));
  g1(12,14)=(-(exp(y(14))*params(21)*exp(y(35))^(-params(22))));
  g1(12,33)=exp(y(33));
  g1(12,35)=(-(exp(y(14))*params(21)*exp(y(35))*getPowerDeriv(exp(y(35)),(-params(22)),1)));
  g1(13,14)=(-T180);
  g1(13,34)=exp(y(34));
  g1(13,36)=(-(exp(y(14))*(1-params(21))*exp(y(36))*getPowerDeriv(exp(y(36)),(-params(22)),1)));
  g1(14,35)=(-(params(21)*exp(y(35))*getPowerDeriv(exp(y(35)),1-params(22),1)));
  g1(14,36)=(-((1-params(21))*exp(y(36))*getPowerDeriv(exp(y(36)),1-params(22),1)));
  g1(15,18)=(-(exp(y(18))*y(23)/(y(23)-1)));
  g1(15,23)=(-(exp(y(18))*(y(23)-1-y(23))/((y(23)-1)*(y(23)-1))));
  g1(15,35)=exp(y(35));
  g1(16,36)=exp(y(36));
  g1(16,37)=(-T197);
  g1(17,29)=(-exp(y(29)));
  g1(17,37)=exp(y(37));
  g1(18,15)=(-(T207*exp(y(26))*exp(y(45))*exp(y(15))*params(30)/exp(y(29))*T391*T394));
  g1(18,24)=(-(T221*exp(y(47))*(-(exp(y(24))*exp(y(12))))/(exp(y(24))*exp(y(24)))*getPowerDeriv(T204,params(15),1)));
  g1(18,26)=(-(T207*T219*T394));
  g1(18,29)=(-(T207*T394*exp(y(26))*T391*(-(exp(y(29))*exp(y(45))*(1+params(30)*(exp(y(15))-1))))/(exp(y(29))*exp(y(29)))));
  g1(18,12)=(-(T221*exp(y(47))*T204*getPowerDeriv(T204,params(15),1)));
  g1(18,42)=exp(y(42));
  g1(18,45)=(-(T207*T394*exp(y(26))*T213*T391));
  g1(18,47)=(-(T207*T221));
  g1(19,35)=(-(exp(y(42))*params(34)*T225*getPowerDeriv(T225,(-params(35)),1)));
  g1(19,42)=(-(exp(y(42))*params(34)*T225^(-params(35))));
  g1(19,43)=exp(y(43));
  g1(19,45)=(-(exp(y(42))*params(34)*getPowerDeriv(T225,(-params(35)),1)*(-(exp(y(45))*exp(y(35))))/(exp(y(45))*exp(y(45)))));
  g1(20,42)=(-T238);
  g1(20,44)=exp(y(44));
  g1(20,45)=(-(exp(y(42))*(1-params(34))*(-(exp(y(45))*exp(y(46))))/(exp(y(45))*exp(y(45)))*getPowerDeriv(T235,(-params(35)),1)));
  g1(20,46)=(-(exp(y(42))*(1-params(34))*T235*getPowerDeriv(T235,(-params(35)),1)));
  g1(21,35)=(-(params(34)*exp(y(35))*getPowerDeriv(exp(y(35)),1-params(35),1)*T550));
  g1(21,45)=exp(y(45));
  g1(21,46)=(-(T550*(1-params(34))*exp(y(46))*getPowerDeriv(exp(y(46)),1-params(35),1)));
  g1(22,37)=(-exp(y(37)));
  g1(22,46)=exp(y(46));
  g1(23,21)=(-T258);
  g1(23,35)=(-(exp(y(21))*params(23)*T253*getPowerDeriv(T253,(-params(24)),1)));
  g1(23,38)=exp(y(38));
  g1(23,40)=(-(exp(y(21))*params(23)*getPowerDeriv(T253,(-params(24)),1)*(-(exp(y(35))*exp(y(40))))/(exp(y(40))*exp(y(40)))));
  g1(24,21)=(-T266);
  g1(24,39)=exp(y(39));
  g1(24,40)=(-(exp(y(21))*(1-params(23))*(-(exp(y(40))*exp(y(41))))/(exp(y(40))*exp(y(40)))*getPowerDeriv(T263,(-params(24)),1)));
  g1(24,41)=(-(exp(y(21))*(1-params(23))*T263*getPowerDeriv(T263,(-params(24)),1)));
  g1(25,35)=(-(params(23)*exp(y(35))*getPowerDeriv(exp(y(35)),1-params(24),1)*T561));
  g1(25,40)=exp(y(40));
  g1(25,41)=(-(T561*(1-params(23))*exp(y(41))*getPowerDeriv(exp(y(41)),1-params(24),1)));
  g1(26,37)=(-T197);
  g1(26,41)=exp(y(41));
  g1(27,6)=(-params(61));
  g1(27,26)=1;
  g1(27,7)=(-params(62));
  g1(27,8)=(-params(63));
  g1(27,62)=(-1);
  g1(28,6)=(-params(64));
  g1(28,7)=(-params(65));
  g1(28,27)=1;
  g1(28,8)=(-params(66));
  g1(28,63)=(-1);
  g1(29,6)=(-params(67));
  g1(29,7)=(-params(68));
  g1(29,8)=(-params(69));
  g1(29,28)=1;
  g1(29,64)=(-1);
  g1(30,4)=(-params(36));
  g1(30,24)=1;
  g1(30,53)=(-1);
  g1(31,3)=(-params(49));
  g1(31,23)=1;
  g1(31,67)=(-1);
  g1(32,5)=(-params(47));
  g1(32,25)=1;
  g1(32,66)=(-1);
  g1(33,11)=(-params(48));
  g1(33,32)=1;
  g1(33,65)=(-1);
  g1(34,13)=(-params(43));
  g1(34,47)=1;
  g1(34,59)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],34,4624);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],34,314432);
end
end
