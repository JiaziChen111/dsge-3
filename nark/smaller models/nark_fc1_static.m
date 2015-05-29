function [residual, g1, g2] = nark_fc1_static(y, x, params)
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

residual = zeros( 19, 1);

%
% Model equations
%

T46 = exp(y(1))-exp(y(1))*params(5)/exp(y(11));
T64 = T46*params(29)/(params(29)-1)*exp(y(4))^params(25);
T92 = exp(y(6))*exp(y(9))*exp(y(12))*(1-params(4))/params(4);
T94 = T92/exp(y(4))/exp(y(11));
T96 = exp(y(3))/(1-params(4));
T99 = (exp(y(9))/params(4))^params(4);
T106 = exp(y(6))^params(4)*exp(y(4))^(1-params(4));
T108 = exp(y(11))^(-params(4));
T147 = exp((-params(14))*(exp(y(16))*exp(y(18))-exp((y(7)))*params(13))-params(33)*(exp(y(14))-exp((y(14)))-(exp(y(2))-exp((y(2))))));
T259 = (-((-(exp(y(11))*exp(y(1))*params(5)))/(exp(y(11))*exp(y(11)))));
lhs =y(11);
rhs =y(11)*params(36)+(1-params(36))*params(18)+x(3);
residual(1)= lhs-rhs;
lhs =y(10);
rhs =(1-params(49))*params(6)+y(10)*params(49)+x(17);
residual(2)= lhs-rhs;
lhs =y(12);
rhs =y(12)*params(47)+x(16);
residual(3)= lhs-rhs;
lhs =y(19);
rhs =y(19)*params(48)+x(15);
residual(4)= lhs-rhs;
lhs =exp(y(11))/T46;
rhs =params(1)/T46*exp(y(2));
residual(5)= lhs-rhs;
lhs =exp(y(3));
rhs =T64;
residual(6)= lhs-rhs;
lhs =1/exp(y(5));
rhs =y(10)/(y(10)-1);
residual(7)= lhs-rhs;
lhs =exp(y(9));
rhs =exp(y(2))-(1-params(2));
residual(8)= lhs-rhs;
lhs =exp(y(6))-(1-params(2))*exp(y(6))/exp(y(11));
rhs =exp(y(8));
residual(9)= lhs-rhs;
lhs =exp(y(3));
rhs =T94;
residual(10)= lhs-rhs;
lhs =exp(y(5));
rhs =T96^(1-params(4))*T99;
residual(11)= lhs-rhs;
lhs =exp(y(7));
rhs =T106*T108;
residual(12)= lhs-rhs;
lhs =exp(y(7));
rhs =exp(y(1))+exp(y(8));
residual(13)= lhs-rhs;
lhs =(exp(y(9))+1-params(2))*exp(y(15))*exp(y(16));
rhs =exp(y(16))*exp(y(17));
residual(14)= lhs-rhs;
lhs =exp(y(17));
rhs =exp(y(14))*T147*exp(y(19));
residual(15)= lhs-rhs;
lhs =exp(y(16))*exp(y(17))*exp(y(18));
rhs =exp(y(15))*exp(y(11))*exp(y(16))*exp(y(18));
residual(16)= lhs-rhs;
lhs =y(13);
rhs =y(13)*params(61)+y(14)*params(62)+y(15)*params(63)+x(12);
residual(17)= lhs-rhs;
lhs =y(14);
rhs =y(13)*params(64)+y(14)*params(65)+y(15)*params(66)+x(13);
residual(18)= lhs-rhs;
lhs =y(15);
rhs =y(13)*params(67)+y(14)*params(68)+y(15)*params(69)+x(14);
residual(19)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(19, 19);

  %
  % Jacobian matrix
  %

  g1(1,11)=1-params(36);
  g1(2,10)=1-params(49);
  g1(3,12)=1-params(47);
  g1(4,19)=1-params(48);
  g1(5,1)=(-(exp(y(11))*T46))/(T46*T46)-exp(y(2))*(-(T46*params(1)))/(T46*T46);
  g1(5,2)=(-(params(1)/T46*exp(y(2))));
  g1(5,11)=(exp(y(11))*T46-exp(y(11))*T259)/(T46*T46)-exp(y(2))*(-(params(1)*T259))/(T46*T46);
  g1(6,1)=(-T64);
  g1(6,3)=exp(y(3));
  g1(6,4)=(-(T46*params(29)/(params(29)-1)*exp(y(4))*getPowerDeriv(exp(y(4)),params(25),1)));
  g1(6,11)=(-(exp(y(4))^params(25)*params(29)/(params(29)-1)*T259));
  g1(7,5)=(-exp(y(5)))/(exp(y(5))*exp(y(5)));
  g1(7,10)=(-((y(10)-1-y(10))/((y(10)-1)*(y(10)-1))));
  g1(8,2)=(-exp(y(2)));
  g1(8,9)=exp(y(9));
  g1(9,6)=exp(y(6))-(1-params(2))*exp(y(6))/exp(y(11));
  g1(9,8)=(-exp(y(8)));
  g1(9,11)=(-((-(exp(y(11))*(1-params(2))*exp(y(6))))/(exp(y(11))*exp(y(11)))));
  g1(10,3)=exp(y(3));
  g1(10,4)=(-((-(exp(y(4))*T92))/(exp(y(4))*exp(y(4)))/exp(y(11))));
  g1(10,6)=(-T94);
  g1(10,9)=(-T94);
  g1(10,11)=(-((-(exp(y(11))*T92/exp(y(4))))/(exp(y(11))*exp(y(11)))));
  g1(10,12)=(-T94);
  g1(11,3)=(-(T99*T96*getPowerDeriv(T96,1-params(4),1)));
  g1(11,5)=exp(y(5));
  g1(11,9)=(-(T96^(1-params(4))*exp(y(9))/params(4)*getPowerDeriv(exp(y(9))/params(4),params(4),1)));
  g1(12,4)=(-(T108*exp(y(6))^params(4)*exp(y(4))*getPowerDeriv(exp(y(4)),1-params(4),1)));
  g1(12,6)=(-(T108*exp(y(4))^(1-params(4))*exp(y(6))*getPowerDeriv(exp(y(6)),params(4),1)));
  g1(12,7)=exp(y(7));
  g1(12,11)=(-(T106*exp(y(11))*getPowerDeriv(exp(y(11)),(-params(4)),1)));
  g1(13,1)=(-exp(y(1)));
  g1(13,7)=exp(y(7));
  g1(13,8)=(-exp(y(8)));
  g1(14,9)=exp(y(16))*exp(y(9))*exp(y(15));
  g1(14,15)=(exp(y(9))+1-params(2))*exp(y(15))*exp(y(16));
  g1(14,16)=(exp(y(9))+1-params(2))*exp(y(15))*exp(y(16))-exp(y(16))*exp(y(17));
  g1(14,17)=(-(exp(y(16))*exp(y(17))));
  g1(15,2)=(-(exp(y(14))*exp(y(19))*T147*(-(params(33)*(-(exp(y(2))-exp((y(2)))))))));
  g1(15,7)=(-(exp(y(14))*exp(y(19))*T147*(-params(14))*(-(exp((y(7)))*params(13)))));
  g1(15,14)=(-(exp(y(14))*T147*exp(y(19))+exp(y(14))*exp(y(19))*T147*(-(params(33)*(exp(y(14))-exp((y(14))))))));
  g1(15,16)=(-(exp(y(14))*exp(y(19))*T147*(-params(14))*exp(y(16))*exp(y(18))));
  g1(15,17)=exp(y(17));
  g1(15,18)=(-(exp(y(14))*exp(y(19))*T147*(-params(14))*exp(y(16))*exp(y(18))));
  g1(15,19)=(-(exp(y(14))*T147*exp(y(19))));
  g1(16,11)=(-(exp(y(15))*exp(y(11))*exp(y(16))*exp(y(18))));
  g1(16,15)=(-(exp(y(15))*exp(y(11))*exp(y(16))*exp(y(18))));
  g1(16,16)=exp(y(16))*exp(y(17))*exp(y(18))-exp(y(15))*exp(y(11))*exp(y(16))*exp(y(18));
  g1(16,17)=exp(y(16))*exp(y(17))*exp(y(18));
  g1(16,18)=exp(y(16))*exp(y(17))*exp(y(18))-exp(y(15))*exp(y(11))*exp(y(16))*exp(y(18));
  g1(17,13)=1-params(61);
  g1(17,14)=(-params(62));
  g1(17,15)=(-params(63));
  g1(18,13)=(-params(64));
  g1(18,14)=1-params(65);
  g1(18,15)=(-params(66));
  g1(19,13)=(-params(67));
  g1(19,14)=(-params(68));
  g1(19,15)=1-params(69);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],19,361);
end
end
