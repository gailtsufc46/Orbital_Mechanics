function d = OrbMvrBiElliptic( rA, rB, rC, mu )
%
% Compute a bi-elliptic maneuver
%
% Inputs: 
%   rA    Initial orbit radius
%   rB    Intermediate (largest) distance at apogee for transfer
%   rC    Final orbit radius
%   mu    Gravitational constant (km^3/s^2)
%
% Outputs:
%   d     Data structure with transfer times and delta-vs for each segment
%

% circular orbit velocities
vA1 = sqrt(mu/rA);
vC4 = sqrt(mu/rC);

% angular momentum of each transfer ellipse
h2 = sqrt(2*mu)*sqrt(rA*rB/(rA+rB));
h3 = sqrt(2*mu)*sqrt(rB*rC/(rB+rC));

% velocity at start / end of first transfer
vA2 = h2/rA;
vB2 = h2/rB;

% velocity at start / end of second transfer
vB3 = h3/rB;
vC3 = h3/rC;

% delta-v's
dVA = abs( vA2 - vA1 );
dVB = abs( vB3 - vB2 );
dVC = abs( vC4 - vC3 );

% total delta-v
dVTot = dVA + dVB + dVC;

% transfer times
e2 = abs(rA-rB)/(rA+rB);
a2 = h2^2/mu/(1-e2^2);
dT2 = OrbPeriod(a2,mu)/2;

e3 = abs(rB-rC)/(rB+rC);
a3 = h3^2/mu/(1-e3^2);
dT3 = OrbPeriod(a3,mu)/2;

dTTot = dT2+dT3;

% gather all variables into data structure "d" for output
w=who;
for i=1:length(w)
  eval(['d.',w{i},'=',w{i},';'])
end
