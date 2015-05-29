function [residual, g1, g2] = nark_fc_static(y, x, params)
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

residual = zeros( 14, 1);

%
% Model equations
%

T40 = exp(y(1))-exp(y(1))*params(5)/exp(y(10));
T61 = T40*params(29)/(params(29)-1)*exp(y(3))^params(25);
T83 = exp(y(5))*exp(y(8))*exp(y(11))*(1-params(4))/params(4);
T85 = T83/exp(y(3))/exp(y(10));
T87 = exp(y(2))/(1-params(4));
T90 = (exp(y(8))/params(4))^params(4);
T97 = exp(y(5))^params(4)*exp(y(3))^(1-params(4));
T99 = exp(y(10))^(-params(4));
T194 = (-((-(exp(y(10))*exp(y(1))*params(5)))/(exp(y(10))*exp(y(10)))));
lhs =y(10);
rhs =y(10)*params(36)+(1-params(36))*params(18)+x(3);
residual(1)= lhs-rhs;
lhs =y(9);
rhs =(1-params(49))*params(6)+y(9)*params(49)+x(17);
residual(2)= lhs-rhs;
lhs =y(11);
rhs =y(11)*params(47)+x(16);
residual(3)= lhs-rhs;
lhs =exp(y(10))/T40;
rhs =params(1)/T40*(exp(y(8))+1-params(2));
residual(4)= lhs-rhs;
lhs =exp(y(2));
rhs =T61;
residual(5)= lhs-rhs;
lhs =1/exp(y(4));
rhs =y(9)/(y(9)-1);
residual(6)= lhs-rhs;
lhs =exp(y(5))-(1-params(2))*exp(y(5))/exp(y(10));
rhs =exp(y(7));
residual(7)= lhs-rhs;
lhs =exp(y(2));
rhs =T85;
residual(8)= lhs-rhs;
lhs =exp(y(4));
rhs =T87^(1-params(4))*T90;
residual(9)= lhs-rhs;
lhs =exp(y(6));
rhs =T97*T99;
residual(10)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(1))+exp(y(7));
residual(11)= lhs-rhs;
lhs =y(12);
rhs =y(12)*params(61)+params(62)*y(13)+params(63)*y(14)+x(12);
residual(12)= lhs-rhs;
lhs =y(13);
rhs =y(12)*params(64)+y(13)*params(65)+y(14)*params(66)+x(13);
residual(13)= lhs-rhs;
lhs =y(14);
rhs =y(12)*params(67)+y(13)*params(68)+y(14)*params(69)+x(14);
residual(14)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(14, 14);

  %
  % Jacobian matrix
  %

  g1(1,10)=1-params(36);
  g1(2,9)=1-params(49);
  g1(3,11)=1-params(47);
  g1(4,1)=(-(exp(y(10))*T40))/(T40*T40)-(exp(y(8))+1-params(2))*(-(T40*params(1)))/(T40*T40);
  g1(4,8)=(-(params(1)/T40*exp(y(8))));
  g1(4,10)=(exp(y(10))*T40-exp(y(10))*T194)/(T40*T40)-(exp(y(8))+1-params(2))*(-(params(1)*T194))/(T40*T40);
  g1(5,1)=(-T61);
  g1(5,2)=exp(y(2));
  g1(5,3)=(-(T40*params(29)/(params(29)-1)*exp(y(3))*getPowerDeriv(exp(y(3)),params(25),1)));
  g1(5,10)=(-(exp(y(3))^params(25)*params(29)/(params(29)-1)*T194));
  g1(6,4)=(-exp(y(4)))/(exp(y(4))*exp(y(4)));
  g1(6,9)=(-((y(9)-1-y(9))/((y(9)-1)*(y(9)-1))));
  g1(7,5)=exp(y(5))-(1-params(2))*exp(y(5))/exp(y(10));
  g1(7,7)=(-exp(y(7)));
  g1(7,10)=(-((-(exp(y(10))*(1-params(2))*exp(y(5))))/(exp(y(10))*exp(y(10)))));
  g1(8,2)=exp(y(2));
  g1(8,3)=(-((-(exp(y(3))*T83))/(exp(y(3))*exp(y(3)))/exp(y(10))));
  g1(8,5)=(-T85);
  g1(8,8)=(-T85);
  g1(8,10)=(-((-(exp(y(10))*T83/exp(y(3))))/(exp(y(10))*exp(y(10)))));
  g1(8,11)=(-T85);
  g1(9,2)=(-(T90*T87*getPowerDeriv(T87,1-params(4),1)));
  g1(9,4)=exp(y(4));
  g1(9,8)=(-(T87^(1-params(4))*exp(y(8))/params(4)*getPowerDeriv(exp(y(8))/params(4),params(4),1)));
  g1(10,3)=(-(T99*exp(y(5))^params(4)*exp(y(3))*getPowerDeriv(exp(y(3)),1-params(4),1)));
  g1(10,5)=(-(T99*exp(y(3))^(1-params(4))*exp(y(5))*getPowerDeriv(exp(y(5)),params(4),1)));
  g1(10,6)=exp(y(6));
  g1(10,10)=(-(T97*exp(y(10))*getPowerDeriv(exp(y(10)),(-params(4)),1)));
  g1(11,1)=(-exp(y(1)));
  g1(11,6)=exp(y(6));
  g1(11,7)=(-exp(y(7)));
  g1(12,12)=1-params(61);
  g1(12,13)=(-params(62));
  g1(12,14)=(-params(63));
  g1(13,12)=(-params(64));
  g1(13,13)=1-params(65);
  g1(13,14)=(-params(66));
  g1(14,12)=(-params(67));
  g1(14,13)=(-params(68));
  g1(14,14)=1-params(69);
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],14,196);
end
end
