% This file is modified by Li Xiangyang @ 2012-5-6 to run in Dynare 4.2.0 or above
% Modfied @2017-7 for illustration.

var c, k lab z;
varexo e;

parameters beta theta delta alpha tau rho s;

beta     = 0.987;
theta     = 0.357;
delta     = 0.012;
alpha     = 0.4;
tau     = 2;
rho     = 0.95;
s       = 0.007;

model; 
% define model-local variables to improve readability;
# mar_c= (c^theta*(1-lab)^(1-theta))^(1-tau);
# mar_c1=(c(+1)^theta*(1-lab(+1))^(1-theta))^(1-tau);

%(1) Euler equation
[name='Euler equation']
mar_c/c=beta*(mar_c1/c(+1))*(1+alpha*exp(z)*k(-1)^(alpha-1)*lab^(1-alpha)-delta);

%(2) wage equation
c=theta/(1-theta)*(1-alpha)*exp(z)*k(-1)^alpha*lab^(-alpha)*(1-lab);

%(3)capital accumulation equation
k=exp(z)*k(-1)^alpha*lab^(1-alpha)-c+(1-delta)*k(-1);

%(4) technology shock
z=rho*z(-1)+e;
end;

initval;
k   = 1;
c   = 1;
lab = 0.3;
z   = 0;
e   = 0;
end;

shocks;
var e=s^2;
end;

steady;
resid;
%if periods not specify, there will be no simulations.
stoch_simul(periods=1000,irf=40,order=1) c k lab z; 

%save the simulated data to file
dynasave('simudata.mat');

figure(1)
subplot(2,2,1)
plot(oo_.endo_simul(1,:),'-b','Linewidth',1)
title('Consumption ')

subplot(2,2,2)
plot(oo_.endo_simul(2,:),'-b','Linewidth',1)
title('Capital Stock ')

subplot(2,2,3)
plot(oo_.endo_simul(3,:),'-b','Linewidth',1)
title('Labor')

subplot(2,2,4)
plot(oo_.endo_simul(4,:),'-b','Linewidth',1)
title('Technology Shock')
