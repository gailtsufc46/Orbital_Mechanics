%% Example 7.1 from the textbook

%% Given

mu = 398600.44;

% spacecraft A
hA = 52059;
eA = 0.025724;
iA = 60;
WA = 40;
wA = 30;
thA = 40;

% spacecraft B
hB = 52362;
eB = 0.0072696;
iB = 50;
WB = 40;
wB = 120;
thB = 40;

%% Step 1 - compute position and velocity from orbital elements

coeA = [ hA^2/mu/(1-eA^2), iA*pi/180, WA*pi/180, wA*pi/180, eA, thA*pi/180 ];
[rA,vA] = RVFromCOE( coeA, mu );

coeB = [ hB^2/mu/(1-eB^2), iB*pi/180, WB*pi/180, wB*pi/180, eB, thB*pi/180 ];
[rB,vB] = RVFromCOE( coeB, mu );

%% Step 2 - compute angular velocity vector "omega"
rAsq = rA'*rA;
omega = cross(rA,vA)/ rAsq;

%% Step 3 - compute angular accleration vector "omega-dot"

omegaDot = -2*rA'*vA/rAsq * omega;

%% step 4 - compute ihat direction
rAmag = sqrt(rAsq);
ihat = rA/rAmag;

%% step 5 - compute khat direction
khat = cross(rA,vA)/hA;

%% step 6 - compute jhat direction
jhat = cross(khat,ihat);

%% step 7 - form the Q matrix

Q = [ihat'; jhat'; khat'];

% check norm of matrix - should be 1 ...
if( abs( norm(Q)-1 ) > 1e-9 )
  error('Norm of Q is not 1.')
end

%% step 8 - compute r_REL in ECI

r_REL_ECI = rB - rA;

%% step 9 - compute v_REL in ECI

v_REL_ECI = vB - vA - cross(omega,r_REL_ECI);

%% step 10 - compute a_REL in ECI
rBmag = sqrt(rB'*rB);
aA = -mu*rA/rAmag^3;
aB = -mu*rB/rBmag^3;
a_REL_ECI = aB - aA - cross(omegaDot,r_REL_ECI) ...
  - 2*cross(omega,v_REL_ECI) ...
  - cross( omega, cross(omega,r_REL_ECI) );

%% step 11 - compute r_REL in LVLH

r_REL_LVLH = Q*r_REL_ECI

%% step 12 - compute v_REL in LVLH

v_REL_LVLH = Q*v_REL_ECI

%% step 13 - compute a_REL in LVLH

a_REL_LVLH = Q*a_REL_ECI

