function d = Flyby( muS, muP, rPH, vPH, vSH0, rPer, side, radP, name, showPlot )

%{
  Compute the turn angle and change in velocity from a flyby maneuver.

  Inputs
    muS   Gravitational constant of sun
    muP   Gravitational constant of planet
    rPH   Position of planet (and spacecraft) in heliocentric frame
    vPH   Velocity of planet in heliocentric frame
    vSH0  Initial spacecraft velocity in heliocentric frame
    v0    Initial heliocentric velocity of spacecraft
    rPer  Periapsis distance at flyby
    side  Which side of the planet ('sunlit' or 'dark')
    radP  Radius of the planet (only used for plotting)
    name  Title to be used on plots
    showPlot - Flag (1 or 0) to show plots or not.

  Outputs
    d     Data structure with fields:

          R: planet distance from sun
        vPH: velocity of the planet in the heliocentric frame
         vP: velocity of the planet in the V-S frame
       vRel: relative velocity of spacecraft w.r.t. planet in V-S frame
       vinf: hyperbolic excess velocity (magnitude of vRel)
          e: eccentricity of the hyperbolic trajectory
      delta: turn angle of the hyperbolic trajectory
          a: semimajor axis of the hyperbolic trajectory
     aimRad: aiming radius
       phi1: heading angle of the incoming velocity vector in V-S frame
       phi2: heading angle of the outgoing velocity vector in V-S frame
      vRad0: radial component of spacecraft's initial total velocity
    vTrans0: transverse component of spacecraft's initial total velocity
      vRadF: radial component of spacecraft's final total velocity
    vTransF: transverse component of spacecraft's initial total velocity
         v0: initial velocity vector of spacecraft in S-W frame
         vF: final velocity vector of spacecraft in S-W frame
         dV: delta-v vector in S-W frame
        v0H: initial velocity vector of spacecraft in heliocentric frame 
        vFH: final velocity vector of spacecraft in heliocentric frame 
        dVH: delta-v vector of spacecraft in heliocentric frame 
         h0: initial angular momentum of spacecraft heliocentric orbit
        th0: initial true anomaly of spacecraft heliocentric orbit
         e0: initial eccentricity of spacecraft heliocentric orbit
         a0: initial semimajor axis of spacecraft heliocentric orbit
         hF: final angular momentum of spacecraft heliocentric orbit
        thF: final true anomaly of spacecraft heliocentric orbit
         eF: final eccentricity of spacecraft heliocentric orbit
         aF: final semimajor axis of spacecraft heliocentric orbit
%}

if( nargin < 1 )

  muS = 132.7e9;
  muP = 42830;

  RE  = 150e6;
  rPH = [1.5*RE;0];
  vPH = [0;sqrt(muS/norm(rPH))];
  
  vSH0 = [3;20];
  
  rPer  = 5000;
  
  side = 'sunlit';
  radP = 3390;
  
  d = Flyby( muS, muP, rPH, vPH, vSH0, rPer, side, radP );
  return
  
end

if( nargin < 9 )
  name = 'Flyby';
end
if( nargin < 10 )
  showPlot = 1;
end

% This is a 2D method. For convenience, allow input vectors with 3
% elements, but remove 3rd element. If not zero throw an error.
if( size(rPH,1)==3 )
  if( abs(rPH(3))>eps )
    error('Non-zero 3rd element of rPH input.')
  end
end
if( size(vPH,1)==3 )
  if( abs(vPH(3))>eps )
    error('Non-zero 3rd element of vPH input.')
  end
end
if( size(vSH0,1)==3 )
  if( abs(vSH0(3))>eps )
    error('Non-zero 3rd element of vSH0 input.')
  end
end
rPH = rPH(1:2);
vPH = vPH(1:2);
vSH0 = vSH0(1:2);

% planetary distance from sun
R = norm(rPH);
d.R = R;

% rotation matrix from heliocentric frame to V-S frame
angMomDir = [0;0;1];
radialDir = [rPH;0]/R;
transverseDir = cross(angMomDir,radialDir);
mat = [transverseDir(1:2)';-radialDir(1:2)'];

% planet velocity in V-S frame
d.vPH = vPH;
d.vP = mat*vPH;

% initial spacecraft velocity in V-S frame
v0 = mat*vSH0;

% relative velocity in V-S coordinates
d.vRel = v0 - d.vP;

% hyperbolic excess velocity: vinf
d.vinf = sqrt( d.vRel(1)^2 + d.vRel(2)^2 );

% part b) ... 2 equations (e, delta)
% turn angle: delta
d.e = 1+rPer*d.vinf^2/muP;
d.delta = 2*asin( 1/d.e );

d.a = rPer/(d.e-1);
d.aimRad = d.a*sqrt(d.e^2-1);

% part c) ... 4 equations (phi1, phi2, vRad, vTrans)

% new radial and transverse vel in heliocentric frame
d.phi1 = atan2(d.vRel(2),d.vRel(1));

% if sunlit and leading (negative transverse relatieve velocity),
% or if dark and trailing (positive transverse relatieve velocity),
% then we add the turn angle
if( (strcmpi(side,'sunlit') && d.vRel(1)<0) || ...
    (strcmpi(side,'dark')   && d.vRel(1)>0 )) 
  d.phi2 = d.phi1 + d.delta;
  turnDir = 1;
else
  % otherwise we subtract the turn angle
  d.phi2 = d.phi1 - d.delta;
  turnDir = -1;
end

d.vRad0 = -v0(2);
d.vTrans0 = v0(1);

d.vRadF   = -d.vP(2) - d.vinf*sin(d.phi2);
d.vTransF = d.vP(1) + d.vinf*cos(d.phi2);
d.v0 = v0;
d.vF = [d.vTransF; -d.vRadF];
d.dV = d.vF-d.v0;

d.v0H = mat'*d.v0;
d.vFH = mat'*d.vF;
d.dVH = d.vFH - d.v0H;

% compute initial theta, a, e
d.h0    = R*d.vTrans0;
ecosth  = d.h0^2 / muS /R - 1;
esinth  = d.vRad0*d.h0/muS;
d.th0   = atan2(esinth,ecosth);
d.e0    = esinth/sin(d.th0);
d.a0    = (d.h0^2/muS)/(1-d.e0^2);


% compute new theta, a, e
d.hF    = R*d.vTransF;
ecosth  = d.hF^2 / muS /R - 1;
esinth  = d.vRadF*d.hF/muS;
d.thF   = atan2(esinth,ecosth);
d.eF    = esinth/sin(d.thF);
d.aF    = (d.hF^2/muS)/(1-d.eF^2);



% hyperbolic trajectory in uv / us coordinates
thinf = acos(-1/d.e);
if( turnDir > 0 )
  thx = linspace(-thinf+.1,thinf-.1);
else
  thx = linspace(thinf-.1,-thinf+.1);
end
rx = d.a*(d.e^2-1)./(1+d.e*cos(thx));
if( turnDir<0 )
  ang = -(d.phi1+pi-thinf);
else
  ang = -(d.phi1+pi+thinf);
end
rotMat = [cos(ang) sin(ang);-sin(ang) cos(ang)];
xy_vs = rotMat*[rx.*cos(thx);rx.*sin(thx)];
dist = max(rx);

if( ~showPlot )
  return
end

figure
plot([0 -cos(d.phi1)]*dist,[0 -sin(d.phi1)]*dist,'b--')
hold on
plot([0 cos(d.phi2)]*dist,[0 sin(d.phi2)]*dist,'r--')
plot(xy_vs(1,:),xy_vs(2,:),'k')
axis equal
grid on
title(name)

thp = linspace(0,2*pi);
xp = radP*cos(thp);
yp = radP*sin(thp);
fill(xp,yp,'b');

%plot(xy_vs(1,2),xy_vs(2,2),'r.','markersize',20)
dxy1 = xy_vs(:,3)-xy_vs(:,2);
quiver(xy_vs(1,2),xy_vs(2,2),dxy1(1),dxy1(2),'color','k','maxheadsize',2)
dxy2 = xy_vs(:,end-1)-xy_vs(:,end-2);
quiver(xy_vs(1,end-2),xy_vs(2,end-2),dxy2(1),dxy2(2),'color','k','maxheadsize',2)

