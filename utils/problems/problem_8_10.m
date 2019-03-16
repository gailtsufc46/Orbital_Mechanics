%% Problem 8.10 from the textbook

mu = 398600.44;

rP = 7000;
vP = 9;

h = rP*vP;
e = vP^2*rP/mu - 1;
a = rP/(1-e);
rA = a*(1+e);
vA = h/rA;

vP2 = 9-.001;
h2 = rP*vP2;
e2 = (vP2^2)*rP/mu - 1;
a2 = rP/(1-e2);
rA2 = a2*(1+e2);
vA2 = h2/rA2;

hf  = @(x) rP*x;
ef  = @(x) (x.^2)*rP/mu-1;
af  = @(x) rP./(1-ef(x));
rAf = @(x) af(x).*(1+ef(x));
vAf = @(x) hf(x) ./ rAf(x);

(-rP*vP^2-2*mu)/rP/vP^2



