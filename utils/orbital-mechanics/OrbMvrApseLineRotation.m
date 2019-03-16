function d = OrbMvrApseLineRotation( rP1, rA1, rP2, rA2, eta, mu, draw )

% Compute an orbit maneuver for a Hohmann Transfer. Choices:
%   Go from A on orbit 1 to B on orbit 2
%   Go from A' on orbit 1 to B' on orbit 2
%
% Inputs: 
%   rP1   Initial orbit perigee
%   rA1   Initial orbit apogee
%   rP2   Final orbit perigee
%   rA2   Final orbit apogee
%   eta   Angle by which to rotate the apse line
%   mu    Gravitational constant (km^3/s^2)
%   draw  (1|0) Draw a plot or not
%
% Outputs:
%   d     Data structure with fields:
%         .dVA  Delta-v at point A
%         .dVB  Delta-v at point B
%         .aT   Semi major axis of transfer orbit
%         .eT   Eccentricity of transfer orbit
%         .dT   Elapsed time for maneuver
%

if( nargin<7 )
  draw = 0;
end

e1 = (rA1-rP1)/(rA1+rP1);
e2 = (rA2-rP2)/(rA2+rP2);

h1 = sqrt(2*mu)*sqrt( rA1*rP1/(rA1+rP1) );
h2 = sqrt(2*mu)*sqrt( rA2*rP2/(rA2+rP2) );

a1 = .5*(rA1+rP1);
a2 = .5*(rA2+rP2);

% a1*(1-e1^2) - h1^2/mu
% a2*(1-e2^2) - h2^2/mu

a = e1*h2^2-e2*h1^2*cos(eta);
b = -e2*h1^2*sin(eta);
c = h1^2-h2^2;

phi = atan(b/a);
theta1 = phi + [1,-1] * acos(c/a*cos(phi));

theta2 = theta1-eta;

r = h1^2/mu./(1+e1*cos(theta1));
vT1 = h1./r;
vR1 = mu/h1*e1*sin(theta1);
gam1 = atan2(vR1,vT1);
v1 = sqrt(vR1.^2 + vT1.^2);

vT2 = h2./r;
vR2 = mu/h2*e2*sin(theta1-eta);
gam2 = atan2(vR2,vT2);
v2 = sqrt(vR2.^2 + vT2.^2);

v1v = [vR1; vT1];
v2v = [vR2; vT2];
dVv = [ norm(v2v(:,1)-v1v(:,1)), norm(v2v(:,2)-v1v(:,2)) ]

dV = sqrt( v1.^2 + v2.^2 - 2*v1.*v2.*cos(gam2-gam1) );

phi = atan2(vR2-vR1,vT2-vT1);

d.r = r;
d.gam1 = gam1;
d.gam2 = gam2;
d.vT1 = vT1;
d.vR1 = vR1;
d.vT2 = vT2;
d.vR2 = vR2;
d.phi = phi;
d.dV = dV;
d.theta1 = theta1;
d.v1 = v1;
d.v2 = v2;


if( ~draw )
  return
end

figure
thx = linspace(0,2*pi);
r1 = a1*(1-e1^2)./(1+e1*cos(thx));
x1 = r1.*cos(thx); 
y1 = r1.*sin(thx);

mat = [cos(eta), -sin(eta); sin(eta), cos(eta)];

r2 = a2*(1-e2^2)./(1+e2*cos(thx));
x2 = r2.*cos(thx);
y2 = r2.*sin(thx);

pos2r = mat*[x2;y2];
x2 = pos2r(1,:);
y2 = pos2r(2,:);

plot(0,0,'k.','markersize',16)
hold on

plot(x1,y1,x2,y2)
axis equal, grid on, zoom on, 
plot(r(1)*cos(theta1(1)),r(1)*sin(theta1(1)),'ro')
plot(r(2)*cos(theta1(2)),r(2)*sin(theta1(2)),'rx')

plot(1.2*[-rA1 rP1], [0 0], 'b--')
plot(1.2*[-rA2 rP2]*cos(eta), 1.2*[-rA2 rP2]*sin(eta), 'r--')









