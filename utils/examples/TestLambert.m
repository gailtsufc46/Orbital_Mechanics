%% Tests for solving Lamberts problem

global mu
mu = 398600.44;


%% test set 1
r1 = [8000;-3000;0];
r2 = [9000;2000;3000];
dT = 15*60;
mu = 398600.44;

% mine
[v1,v2]         = LambertSolver(r1,r2,dT,mu,'pro');
[a,i,W,w,e,th]  = OrbitalElementsFromRV(r1,v1,mu);

% textbook
[v1t,v2t] = lambert(r1,r2,dT);
elem = coe_from_sv(r1,v1,mu);
%[h e RA incl w TA a]
at = elem(7);
it = elem(4);
Wt = elem(3);
wt = elem(5);
et = elem(2); 
tht = elem(6);

% retrograde
% [v1b,v2b] = lambert(r1,r2,dT,'retro');
% [ab,ib,Wb,wb,eb,thb]  = OrbitalElementsFromRV(r1,v1b,mu)

%% test set 2
r1 = [0;1000;12000];
r2 = [500;8000;-500];
dT = 105*60;
mu = 398600.44;

% mine
[v1,v2]         = LambertSolver(r1,r2,dT,mu,'pro');
[a,i,W,w,e,th]  = OrbitalElementsFromRV(r1,v1,mu);

% textbook
[v1t,v2t] = lambert(r1,r2,dT);
elem = coe_from_sv(r1,v1t,mu);
%[h e RA incl w TA a]
at = elem(7);
it = elem(4);
Wt = elem(3);
wt = elem(5);
et = elem(2); 
tht = elem(6);

% retrograde
% [v1b,v2b] = lambert(r1,r2,dT,'retro');
% [ab,ib,Wb,wb,eb,thb]  = OrbitalElementsFromRV(r1,v1b,mu)

%% test set 3

mu = 398600.44;

a = 18000; 
inc = 45*pi/180; 
W = 45*pi/180;
w = 90*pi/180;
e = 0.5;
th = 60*pi/180;
coe = [a, inc, W, w, e, th];
[r1,v1] = RVFromCOE(a, inc, W, w, e, th, mu );

dT = 5*60*60;
[r2,v2] = RVAtTFromR0V0(r1,v1,dT,mu);

[v1s,v2s] = LambertSolver(r1,r2,dT,mu,'pro');
[as,incs,Ws,ws,es,ths]  = OrbitalElementsFromRV(r1,v1s,mu);

%% test set 4 - hyperbolic!

mu = 398600.44;

a = -50000; 
inc = 60*pi/180; 
W = 15*pi/180;
w = 10*pi/180;
e = 1.5;
th = 0;
coe = [a, inc, W, w, e, th];
[r1,v1] = RVFromCOE(a, inc, W, w, e, th, mu );

dT = 60*60;
[r2,v2] = RVAtTFromR0V0(r1,v1,dT,mu);

[v1s,v2s] = LambertSolver(r1,r2,dT,mu,'pro');
[as,incs,Ws,ws,es,ths]  = OrbitalElementsFromRV(r1,v1s,mu);



