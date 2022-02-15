% This file is modified by Li Xiangyang @ 2013-8-30 to run in Dynare 4.2.0 or above
% Modified at 2017-7@Shanghai for illustration.

var y c k inv lab z;
varexo e;

parameters bet the del alp tau rho s;

bet     = 0.987;
the     = 0.357;
del     = 0.012;
alp     = 0.4;
tau     = 2;
rho     = 0.95;
s       = 0.007;

model; 
% define model-local variables to improve readability;
# mar_c = (exp(c*the)*(1-exp(lab))^(1-the))^(1-tau);
# mar_c1= (exp(c(+1)*the)*(1-exp(lab(+1)))^(1-the))^(1-tau);

%(1) Euler equation
mar_c/exp(c)=bet*(mar_c1/exp(c(+1)))*(1+alp*exp(z)*exp(k(-1))^(alp-1)*exp(lab)^(1-alp)-del);

%(2) wage equation
exp(c)=the/(1-the)*(1-alp)*exp(z)*exp(k(-1))^alp*exp(lab)^(-alp)*(1-exp(lab));

%(3) capital accumulation equation
exp(k)=exp(inv)+(1-del)*exp(k(-1));

%(4) the production technology
exp(y)=exp(z)*exp(k(-1))^alp*exp(lab)^(1-alp);

%(5) the resource constraint
exp(y)=exp(c)+exp(inv);

%(6)the technology shock
z=rho*z(-1)+e;
end;

initval;
k   = log(29.71828);
c   = log(1.4);
lab = log(0.3);
z   = log(1);
e   = log(1);
end;

steady;
check;

shocks;
var e=s^2;
end;



%if periods not specify, there will be no simulations.
stoch_simul(periods=1000,irf=40,order=1) c k lab z; 

%save the simulated data to file
dynasave('simudata_exp.mat');
