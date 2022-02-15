function [residual, g1, g2, g3] = GrowthApproximate_dynamic(y, x, params, steady_state, it_)
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

residual = zeros(4, 1);
mar_c__ = (y(3)^params(2)*(1-y(5))^(1-params(2)))^(1-params(5));
mar_c1__ = (y(7)^params(2)*(1-y(8))^(1-params(2)))^(1-params(5));
T31 = params(1)*mar_c1__/y(7);
T41 = y(5)^(1-params(4));
T45 = 1+params(4)*exp(y(6))*y(1)^(params(4)-1)*T41-params(3);
T52 = exp(y(6))*(1-params(4))*params(2)/(1-params(2))*y(1)^params(4);
T55 = T52*y(5)^(-params(4));
T74 = getPowerDeriv(y(3)^params(2)*(1-y(5))^(1-params(2)),1-params(5),1);
T82 = getPowerDeriv(y(7)^params(2)*(1-y(8))^(1-params(2)),1-params(5),1);
lhs =mar_c__/y(3);
rhs =T31*T45;
residual(1)= lhs-rhs;
lhs =y(3);
rhs =(1-y(5))*T55;
residual(2)= lhs-rhs;
lhs =y(4);
rhs =T41*exp(y(6))*y(1)^params(4)-y(3)+y(1)*(1-params(3));
residual(3)= lhs-rhs;
lhs =y(6);
rhs =params(6)*y(2)+x(it_, 1);
residual(4)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(4, 9);

  %
  % Jacobian matrix
  %

  g1(1,3)=(y(3)*(1-y(5))^(1-params(2))*getPowerDeriv(y(3),params(2),1)*T74-mar_c__)/(y(3)*y(3));
  g1(1,7)=(-(T45*params(1)*(y(7)*(1-y(8))^(1-params(2))*getPowerDeriv(y(7),params(2),1)*T82-mar_c1__)/(y(7)*y(7))));
  g1(1,1)=(-(T31*T41*params(4)*exp(y(6))*getPowerDeriv(y(1),params(4)-1,1)));
  g1(1,5)=T74*y(3)^params(2)*(-(getPowerDeriv(1-y(5),1-params(2),1)))/y(3)-T31*params(4)*exp(y(6))*y(1)^(params(4)-1)*getPowerDeriv(y(5),1-params(4),1);
  g1(1,8)=(-(T45*params(1)*T82*y(7)^params(2)*(-(getPowerDeriv(1-y(8),1-params(2),1)))/y(7)));
  g1(1,6)=(-(T31*params(4)*exp(y(6))*y(1)^(params(4)-1)*T41));
  g1(2,3)=1;
  g1(2,1)=(-((1-y(5))*y(5)^(-params(4))*exp(y(6))*(1-params(4))*params(2)/(1-params(2))*getPowerDeriv(y(1),params(4),1)));
  g1(2,5)=(-((-T55)+(1-y(5))*T52*getPowerDeriv(y(5),(-params(4)),1)));
  g1(2,6)=(-((1-y(5))*T55));
  g1(3,3)=1;
  g1(3,1)=(-(1-params(3)+T41*exp(y(6))*getPowerDeriv(y(1),params(4),1)));
  g1(3,4)=1;
  g1(3,5)=(-(exp(y(6))*y(1)^params(4)*getPowerDeriv(y(5),1-params(4),1)));
  g1(3,6)=(-(T41*exp(y(6))*y(1)^params(4)));
  g1(4,2)=(-params(6));
  g1(4,6)=1;
  g1(4,9)=(-1);
end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],4,81);
end
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],4,729);
end
end
