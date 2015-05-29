function [residual, g1, g2] = nark_fc_fb_xm_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                     columns: variables in declaration order
%                                                     rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 34, 1);

%
% Model equations
%

T15 = exp(y(1))-exp(y(1))*params(5)/exp(y(11));
T33 = T15*params(29)/(params(29)-1)*exp(y(4))^params(25);
T56 = exp(y(6))*exp(y(9))*exp(y(12))*(1-params(4))/params(4);
T58 = T56/exp(y(4))/exp(y(11));
T62 = exp(y(3))/(1-params(4));
T65 = (exp(y(9))/params(4))^params(4);
T72 = exp(y(6))^params(4)*exp(y(4))^(1-params(4));
T74 = exp(y(11))^(-params(4));
T120 = exp((-params(14))*(exp(y(16))*exp(y(18))-exp((y(7)))*params(13))-params(33)*(exp(y(14))-exp((y(14)))-(exp(y(2))-exp((y(2))))));
T132 = exp(y(16))*exp(y(17))*exp(y(18))/exp(y(11));
T162 = exp(y(1))*(1-params(21))*exp(y(23))^(-params(22));
T179 = exp(y(24))*params(7)/(params(7)-1);
T184 = exp(y(29))/exp(y(11));
T187 = exp(y(34))*T184^params(15);
T193 = exp(y(32))*(1+params(30)*(exp(y(2))-1))/exp(y(16));
T199 = T193^(-params(16))*exp(y(13));
T201 = T199^(1-params(15));
T205 = exp(y(22))/exp(y(32));
T215 = exp(y(33))/exp(y(32));
T218 = exp(y(29))*(1-params(34))*T215^(-params(35));
T225 = params(34)*exp(y(22))^(1-params(35))+(1-params(34))*exp(y(33))^(1-params(35));
T233 = exp(y(22))/exp(y(27));
T238 = exp(y(8))*params(23)*T233^(-params(24));
T243 = exp(y(28))/exp(y(27));
T246 = exp(y(8))*(1-params(23))*T243^(-params(24));
T253 = params(23)*exp(y(22))^(1-params(24))+(1-params(23))*exp(y(28))^(1-params(24));
T348 = getPowerDeriv(T193,(-params(16)),1);
T351 = getPowerDeriv(T199,1-params(15),1);
T405 = (-((-(exp(y(11))*exp(y(1))*params(5)))/(exp(y(11))*exp(y(11)))));
T505 = getPowerDeriv(T225,1/(1-params(35)),1);
T516 = getPowerDeriv(T253,1/(1-params(24)),1);
lhs =exp(y(11))/T15;
rhs =params(1)/T15*exp(y(2));
residual(1)= lhs-rhs;
lhs =exp(y(3));
rhs =T33;
residual(2)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(y(2))-(1-params(2));
residual(3)= lhs-rhs;
lhs =exp(y(6))-(1-params(2))*exp(y(6))/exp(y(11));
rhs =exp(y(8));
residual(4)= lhs-rhs;
lhs =exp(y(3));
rhs =T58;
residual(5)= lhs-rhs;
lhs =exp(y(5));
rhs =T62^(1-params(4))*T65;
residual(6)= lhs-rhs;
lhs =exp(y(7));
rhs =T72*T74;
residual(7)= lhs-rhs;
lhs =exp(y(7));
rhs =exp(y(20))+exp(y(25))+exp(y(30));
residual(8)= lhs-rhs;
lhs =(exp(y(9))+1-params(2))*exp(y(15))*exp(y(16));
rhs =exp(y(16))*exp(y(17));
residual(9)= lhs-rhs;
lhs =exp(y(17));
rhs =exp(y(14))*T120*exp(y(19));
residual(10)= lhs-rhs;
lhs =exp(y(32))*exp(y(29))+T132/exp(y(15));
rhs =exp(y(16))*exp(y(18))+exp(y(24))*(exp(y(21))+exp(y(31))+exp(y(26)));
residual(11)= lhs-rhs;
lhs =exp(y(20));
rhs =exp(y(1))*params(21)*exp(y(22))^(-params(22));
residual(12)= lhs-rhs;
lhs =exp(y(21));
rhs =T162;
residual(13)= lhs-rhs;
lhs =1;
rhs =params(21)*exp(y(22))^(1-params(22))+(1-params(21))*exp(y(23))^(1-params(22));
residual(14)= lhs-rhs;
lhs =exp(y(22));
rhs =exp(y(5))*y(10)/(y(10)-1);
residual(15)= lhs-rhs;
lhs =exp(y(23));
rhs =T179;
residual(16)= lhs-rhs;
lhs =exp(y(24));
rhs =exp(y(16));
residual(17)= lhs-rhs;
lhs =exp(y(29));
rhs =T187*T201;
residual(18)= lhs-rhs;
lhs =exp(y(30));
rhs =exp(y(29))*params(34)*T205^(-params(35));
residual(19)= lhs-rhs;
lhs =exp(y(31));
rhs =T218;
residual(20)= lhs-rhs;
lhs =exp(y(32));
rhs =T225^(1/(1-params(35)));
residual(21)= lhs-rhs;
lhs =exp(y(33));
rhs =exp(y(24));
residual(22)= lhs-rhs;
lhs =exp(y(25));
rhs =T238;
residual(23)= lhs-rhs;
lhs =exp(y(26));
rhs =T246;
residual(24)= lhs-rhs;
lhs =exp(y(27));
rhs =T253^(1/(1-params(24)));
residual(25)= lhs-rhs;
lhs =exp(y(28));
rhs =T179;
residual(26)= lhs-rhs;
lhs =y(13);
rhs =y(13)*params(61)+y(14)*params(62)+y(15)*params(63)+x(12);
residual(27)= lhs-rhs;
lhs =y(14);
rhs =y(13)*params(64)+y(14)*params(65)+y(15)*params(66)+x(13);
residual(28)= lhs-rhs;
lhs =y(15);
rhs =y(13)*params(67)+y(14)*params(68)+y(15)*params(69)+x(14);
residual(29)= lhs-rhs;
lhs =y(11);
rhs =y(11)*params(36)+(1-params(36))*params(18)+x(3);
residual(30)= lhs-rhs;
lhs =y(10);
rhs =(1-params(49))*params(6)+y(10)*params(49)+x(17);
residual(31)= lhs-rhs;
lhs =y(12);
rhs =y(12)*params(47)+x(16);
residual(32)= lhs-rhs;
lhs =y(19);
rhs =y(19)*params(48)+x(15);
residual(33)= lhs-rhs;
lhs =y(34);
rhs =y(34)*params(43)+x(9);
residual(34)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(34, 34);

  %
  % Jacobian matrix
  %

  g1(1,1)=(-(exp(y(11))*T15))/(T15*T15)-exp(y(2))*(-(T15*params(1)))/(T15*T15);
  g1(1,2)=(-(params(1)/T15*exp(y(2))));
  g1(1,11)=(exp(y(11))*T15-exp(y(11))*T405)/(T15*T15)-exp(y(2))*(-(params(1)*T405))/(T15*T15);
  g1(2,1)=(-T33);
  g1(2,3)=exp(y(3));
  g1(2,4)=(-(T15*params(29)/(params(29)-1)*exp(y(4))*getPowerDeriv(exp(y(4)),params(25),1)));
  g1(2,11)=(-(exp(y(4))^params(25)*params(29)/(params(29)-1)*T405));
  g1(3,2)=(-exp(y(2)));
  g1(3,9)=exp(y(9));
  g1(4,6)=exp(y(6))-(1-params(2))*exp(y(6))/exp(y(11));
  g1(4,8)=(-exp(y(8)));
  g1(4,11)=(-((-(exp(y(11))*(1-params(2))*exp(y(6))))/(exp(y(11))*exp(y(11)))));
  g1(5,3)=exp(y(3));
  g1(5,4)=(-((-(exp(y(4))*T56))/(exp(y(4))*exp(y(4)))/exp(y(11))));
  g1(5,6)=(-T58);
  g1(5,9)=(-T58);
  g1(5,11)=(-((-(exp(y(11))*T56/exp(y(4))))/(exp(y(11))*exp(y(11)))));
  g1(5,12)=(-T58);
  g1(6,3)=(-(T65*T62*getPowerDeriv(T62,1-params(4),1)));
  g1(6,5)=exp(y(5));
  g1(6,9)=(-(T62^(1-params(4))*exp(y(9))/params(4)*getPowerDeriv(exp(y(9))/params(4),params(4),1)));
  g1(7,4)=(-(T74*exp(y(6))^params(4)*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(4),1)));
  g1(7,6)=(-(T74*exp(y(4))^(1-params(4))*exp(y(6))*getPowerDeriv(exp(y(6)),params(4),1)));
  g1(7,7)=exp(y(7));
  g1(7,11)=(-(T72*exp(y(11))*getPowerDeriv(exp(y(11)),(-params(4)),1)));
  g1(8,7)=exp(y(7));
  g1(8,20)=(-exp(y(20)));
  g1(8,25)=(-exp(y(25)));
  g1(8,30)=(-exp(y(30)));
  g1(9,9)=exp(y(16))*exp(y(9))*exp(y(15));
  g1(9,15)=(exp(y(9))+1-params(2))*exp(y(15))*exp(y(16));
  g1(9,16)=(exp(y(9))+1-params(2))*exp(y(15))*exp(y(16))-exp(y(16))*exp(y(17));
  g1(9,17)=(-(exp(y(16))*exp(y(17))));
  g1(10,2)=(-(exp(y(14))*exp(y(19))*T120*(-(params(33)*(-(exp(y(2))-exp((y(2)))))))));
  g1(10,7)=(-(exp(y(14))*exp(y(19))*T120*(-params(14))*(-(exp((y(7)))*params(13)))));
  g1(10,14)=(-(exp(y(14))*T120*exp(y(19))+exp(y(14))*exp(y(19))*T120*(-(params(33)*(exp(y(14))-exp((y(14))))))));
  g1(10,16)=(-(exp(y(14))*exp(y(19))*T120*(-params(14))*exp(y(16))*exp(y(18))));
  g1(10,17)=exp(y(17));
  g1(10,18)=(-(exp(y(14))*exp(y(19))*T120*(-params(14))*exp(y(16))*exp(y(18))));
  g1(10,19)=(-(exp(y(14))*T120*exp(y(19))));
  g1(11,11)=(-(exp(y(11))*exp(y(16))*exp(y(17))*exp(y(18))))/(exp(y(11))*exp(y(11)))/exp(y(15));
  g1(11,15)=(-(exp(y(15))*T132))/(exp(y(15))*exp(y(15)));
  g1(11,16)=T132/exp(y(15))-exp(y(16))*exp(y(18));
  g1(11,17)=T132/exp(y(15));
  g1(11,18)=T132/exp(y(15))-exp(y(16))*exp(y(18));
  g1(11,21)=(-(exp(y(24))*exp(y(21))));
  g1(11,24)=(-(exp(y(24))*(exp(y(21))+exp(y(31))+exp(y(26)))));
  g1(11,26)=(-(exp(y(24))*exp(y(26))));
  g1(11,29)=exp(y(32))*exp(y(29));
  g1(11,31)=(-(exp(y(24))*exp(y(31))));
  g1(11,32)=exp(y(32))*exp(y(29));
  g1(12,1)=(-(exp(y(1))*params(21)*exp(y(22))^(-params(22))));
  g1(12,20)=exp(y(20));
  g1(12,22)=(-(exp(y(1))*params(21)*exp(y(22))*getPowerDeriv(exp(y(22)),(-params(22)),1)));
  g1(13,1)=(-T162);
  g1(13,21)=exp(y(21));
  g1(13,23)=(-(exp(y(1))*(1-params(21))*exp(y(23))*getPowerDeriv(exp(y(23)),(-params(22)),1)));
  g1(14,22)=(-(params(21)*exp(y(22))*getPowerDeriv(exp(y(22)),1-params(22),1)));
  g1(14,23)=(-((1-params(21))*exp(y(23))*getPowerDeriv(exp(y(23)),1-params(22),1)));
  g1(15,5)=(-(exp(y(5))*y(10)/(y(10)-1)));
  g1(15,10)=(-(exp(y(5))*(y(10)-1-y(10))/((y(10)-1)*(y(10)-1))));
  g1(15,22)=exp(y(22));
  g1(16,23)=exp(y(23));
  g1(16,24)=(-T179);
  g1(17,16)=(-exp(y(16)));
  g1(17,24)=exp(y(24));
  g1(18,2)=(-(T187*exp(y(13))*exp(y(32))*exp(y(2))*params(30)/exp(y(16))*T348*T351));
  g1(18,11)=(-(T201*exp(y(34))*(-(exp(y(11))*exp(y(29))))/(exp(y(11))*exp(y(11)))*getPowerDeriv(T184,params(15),1)));
  g1(18,13)=(-(T187*T199*T351));
  g1(18,16)=(-(T187*T351*exp(y(13))*T348*(-(exp(y(16))*exp(y(32))*(1+params(30)*(exp(y(2))-1))))/(exp(y(16))*exp(y(16)))));
  g1(18,29)=exp(y(29))-T201*exp(y(34))*T184*getPowerDeriv(T184,params(15),1);
  g1(18,32)=(-(T187*T351*exp(y(13))*T193*T348));
  g1(18,34)=(-(T187*T201));
  g1(19,22)=(-(exp(y(29))*params(34)*T205*getPowerDeriv(T205,(-params(35)),1)));
  g1(19,29)=(-(exp(y(29))*params(34)*T205^(-params(35))));
  g1(19,30)=exp(y(30));
  g1(19,32)=(-(exp(y(29))*params(34)*getPowerDeriv(T205,(-params(35)),1)*(-(exp(y(32))*exp(y(22))))/(exp(y(32))*exp(y(32)))));
  g1(20,29)=(-T218);
  g1(20,31)=exp(y(31));
  g1(20,32)=(-(exp(y(29))*(1-params(34))*(-(exp(y(32))*exp(y(33))))/(exp(y(32))*exp(y(32)))*getPowerDeriv(T215,(-params(35)),1)));
  g1(20,33)=(-(exp(y(29))*(1-params(34))*T215*getPowerDeriv(T215,(-params(35)),1)));
  g1(21,22)=(-(params(34)*exp(y(22))*getPowerDeriv(exp(y(22)),1-params(35),1)*T505));
  g1(21,32)=exp(y(32));
  g1(21,33)=(-(T505*(1-params(34))*exp(y(33))*getPowerDeriv(exp(y(33)),1-params(35),1)));
  g1(22,24)=(-exp(y(24)));
  g1(22,33)=exp(y(33));
  g1(23,8)=(-T238);
  g1(23,22)=(-(exp(y(8))*params(23)*T233*getPowerDeriv(T233,(-params(24)),1)));
  g1(23,25)=exp(y(25));
  g1(23,27)=(-(exp(y(8))*params(23)*getPowerDeriv(T233,(-params(24)),1)*(-(exp(y(22))*exp(y(27))))/(exp(y(27))*exp(y(27)))));
  g1(24,8)=(-T246);
  g1(24,26)=exp(y(26));
  g1(24,27)=(-(exp(y(8))*(1-params(23))*(-(exp(y(27))*exp(y(28))))/(exp(y(27))*exp(y(27)))*getPowerDeriv(T243,(-params(24)),1)));
  g1(24,28)=(-(exp(y(8))*(1-params(23))*T243*getPowerDeriv(T243,(-params(24)),1)));
  g1(25,22)=(-(params(23)*exp(y(22))*getPowerDeriv(exp(y(22)),1-params(24),1)*T516));
  g1(25,27)=exp(y(27));
  g1(25,28)=(-(T516*(1-params(23))*exp(y(28))*getPowerDeriv(exp(y(28)),1-params(24),1)));
  g1(26,24)=(-T179);
  g1(26,28)=exp(y(28));
  g1(27,13)=1-params(61);
  g1(27,14)=(-params(62));
  g1(27,15)=(-params(63));
  g1(28,13)=(-params(64));
  g1(28,14)=1-params(65);
  g1(28,15)=(-params(66));
  g1(29,13)=(-params(67));
  g1(29,14)=(-params(68));
  g1(29,15)=1-params(69);
  g1(30,11)=1-params(36);
  g1(31,10)=1-params(49);
  g1(32,12)=1-params(47);
  g1(33,19)=1-params(48);
  g1(34,34)=1-params(43);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],34,1156);
end
end
