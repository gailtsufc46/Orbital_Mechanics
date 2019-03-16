function d = OrbMvrHohmann( rA1, rB2, rA1x, rB2x, mu )
%
% Compute an orbit maneuver for a Hohmann Transfer. Choices:
%   Go from A on orbit 1 to B on orbit 2
%   Go from A' on orbit 1 to B' on orbit 2
%
% Inputs: 
%   rA1   Radius at point A on initial orbit
%   rA1x  Radius at point A' on initial orbit (180 deg from A)
%   rB2   Radius at point B on final orbit
%   rB2x  Radius at point B' on final orbit (180 deg from B)
%   mu    Gravitational constant (km^3/s^2)
%
% Outputs:
%   d     Data structure with fields:
%         .dVA  Delta-v at point A
%         .dVB  Delta-v at point B
%         .aT   Semi major axis of transfer orbit
%         .eT   Eccentricity of transfer orbit
%         .dT   Elapsed time for maneuver
%

% angular momentum of initial and final orbits
h1  = sqrt(2*mu)*sqrt( rA1*rA1x /(rA1+rA1x ) );
h2  = sqrt(2*mu)*sqrt( rB2*rB2x /(rB2+rB2x ) );
hT  = sqrt(2*mu)*sqrt( rA1*rB2  /(rA1+rB2  ) );
hTx = sqrt(2*mu)*sqrt( rA1x*rB2x/(rA1x+rB2x) );

% velocities for first option, A to B
vA1 = h1/rA1;
vB2 = h2/rB2;
vAT  = hT/rA1;
vBT  = hT/rB2;
dVA  = abs(vAT-vA1);
dVB  = abs(vB2-vBT);

% velocities for second option, A' to B'
vA1x = h1/rA1x;
vB2x = h2/rB2x;
vATx = hTx/rA1x;
vBTx = hTx/rB2x;
dVAx = abs(vATx-vA1x);
dVBx = abs(vB2x-vBTx);

dVTot = dVA+dVB;
dVTotx = dVAx+dVBx;

e = abs(rA1-rB2)/(rA1+rB2);
a = hT^2/mu/(1-e^2);
T = OrbPeriod(a,mu);
dT = T/2;

ex = abs(rA1x-rB2x)/(rA1x+rB2x);
ax = hTx^2/mu/(1-ex^2);
Tx = OrbPeriod(a,mu);
dTx = Tx/2;

% gather all variables into data structure "d" for output
w=who;
for i=1:length(w)
  eval(['d.',w{i},'=',w{i},';'])
end

