%% Problem 6.34 from the textbook
% Big orbit change "chase maneuver" with Lambert's problem solution

%% Given
mu = 398600.44;

r1A = [10000; 0; 0];
r2B = [0; 16000; 0];
tOF = 3600;
e2 = 0.5;

%% Initial velocities
v1A = [0;sqrt(mu/r1A(1));0];

h2  = sqrt(r2B(2)*mu*(1+e2));
v2B = [0;0;h2/r2B(2)];

%% Compute position on orbit 2 at point C

% Can use this method  
r2C = RVAtTFromR0V0(r2B,v2B,tOF,mu);

% or another way ...
a2 = r2B(2)/(1-e2);
th2C = TrueAnomFromTime(tOF,a2,e2,mu);
r2C_mag = a2*(1-e2^2)/(1+e2*cos(th2C));
r2C_2 = r2C_mag * [0; cos(th2C); sin(th2C)];

%% Solve for initial and final velocities using Lambert, at A, C
[v3A,v3C] = LambertSolver( r1A, r2C, tOF, mu, 'pro' );

dVA = norm(v3A-v1A)
