function [f,g,fdot,gdot] = LagrangeCoeff( r0, v0, mu, dTh )

% Compute the Lagrange coefficients. This method implements equation 2.158
% (a,b,c,d) from the Curtis textbook.

if( abs(sin(dTh)) < 1e-8 )
  f = 1;
  g = 0;
  fdot = 0;
  gdot = 1;
  return
end

% compute specific angular momentum
hvec = cross(r0,v0);
h = sqrt(hvec'*hvec);

% terms that appear many times:
c = cos(dTh);
s = sin(dTh);
mu_h = mu/h;
mu_h2 = mu/(h^2);

% compute magnitude of r0
r0mag = sqrt(r0'*r0);

% radial component of initial velocity
vr0 = v0'*r0/r0mag;

% compute "r" magnitude (Eqn. 2.152)
rmag = h*h/mu ./ ( 1 + (h*h/mu/r0mag-1)*c - h/mu*vr0*s );

% Lagrange coefficients
f = 1-mu_h2*rmag.*(1-c);
g = r0mag/h*rmag.*s;
fdot = mu_h * (1-c)./s.*(mu_h2*(1-c)-1/r0mag-1./rmag);
gdot = 1-mu_h2*r0mag.*(1-c);

