%% Execute a flyby maneuver at a planet after a Hohmann transfer 
%
% First we have a Hohmann transfer from planet A to planet B.
% However, we do *NOT* execute the 2nd delta-v of the Hohmann...
% Instead, once we enter the SOI at planet B, we do a flyby maneuver. 
% Note that the flyby maneuver is really just a hyperbolic orbit, where the
% gravity of the planet changes the velocity of the spacecraft.
%
% Assume that: 
%   * we are given the aiming radius of the hyperbolic orbit.
%   * all planets follow a circular orbit with zero inclination. 
%
% We want to find the following information about the spacecraft AFTER the 
% flyby maneuver: 
%   * Total velocity vector in the heliocentric frame
%   * semi-major axis (of its heliocentric orbit)
%   * eccentricity (of its heliocentric orbit)
%
% Strategy:
%   0. Compute the total heliocentric velocity of the planet and the 
%       spacecraft when it is at the end of the Hohmann transfer maneuver.
%   
%   1. Compute the initial relative velocity of the spacecraft with respect
%       to the planet, when the spacecraft is entering the SOI: "vInfVector1"
% 
%   2. Compute the turn angle that corresponds to this hyperbolic orbit: 
%       "delta". This is the angle between the two asymptotes.
%
%   3. Compute the NEW relative velocity of the spacecraft with respect
%       to the planet, when the spacecraft is exiting the SOI: "vInfVector2"
%
%   4. Compute the NEW absolute velocity of the spacecraft and the
%       corresponding ortial elements (SMA, eccentricity) post-flyby.
%
% 

%% Parameters 
planetA = 'earth';
planetB = 'jupiter';
%aimingRadius = 3502082.3; % leads to 200,000 km altitude for Earth-Jupiter
aimingRadius = 3200000; 
side = 'sunlit';


%% Planet data and constants
muEarth = 398600.44;
muS = 1.32712440018*1e11; % sun
au = 149597871; 
planets  = {'mercury','venus','earth','mars','jupiter','saturn','neptune','uranus'};
rPlanets = [0.3871 0.72333 0.99918 1.5236 5.2017 9.5502 19.158 29.981]; % AU
radPlanets = [2440 6051.8 6378.1 3389.9 71492 60268 25559 24766]; % km
mPlanets = [0.0553 0.815 1 0.107 318 95.2 14.5 17.2]; % Earth mass
kA = find(strncmp(planetA,planets,10));
kB = find(strncmp(planetB,planets,10));
rA = rPlanets(kA)*au; 
rB = rPlanets(kB)*au;
muP = muEarth*mPlanets(kB); % mu value of planet B (km^3/sec^2)
radP = radPlanets(kB);      % radius of planet B (km)

%% Hohmann transfer...

% angular momentum of the transfer orbit
hT = sqrt(2*muS)*sqrt(rA*rB/(rA+rB));

% velocity of the transfer orbit when it reaches point B
vTB = hT/rB;

% velocity of planet B 
vB = sqrt(muS/rB);

%% Compute the velocity relative to the planet when we enter the SOI

% We will use a "V-S" coordinate system, where:
%   V points along the planet velocity vector, and
%   S points toward the sun.
velVectorOfPlanet = [vB; 0];

% Because this is a Hohmann transfer... the spacecraft velocity is aligned
% in the same direction as the planet.
velVectorOfSpacecraft = [vTB; 0];

% The relative velocity of the spacecraft with respect to the planet is 
% just the vector difference. This is the "v-infinity" vector:
vInfVector1 = velVectorOfSpacecraft - velVectorOfPlanet;
phi1 = atan2(vInfVector1(2),vInfVector1(1));

%% Compute the radius at periapse and the turn angle
vinf = norm(vInfVector1);
% solve Eq. 8.57 for rp. It is a quadratic equation, so it has two roots.
rpSol = roots( [1, 2*muP/vinf, -aimingRadius^2] );
% If this problem is feasible, it will have two real roots (one pos, one neg). 
% Pick the positive one!
rp = max(rpSol);
fprintf(1,'Flyby altitude: %5.1f km\n',rp-radP)
% now we can compute the eccentricity and the turn angle.
e = 1 + rp*vinf^2/muP;
delta = 2*asin(1/e);

%% Compute the velocity relative to the planet when we exit the SOI
% It has the same magnitude, but it has rotated by "delta".
% Think about the different geometries that are possible:
%   * SUNLIT and SLOW -- spacecraft flies by the planet on the sunlit side,
%                         and the spacecraft is slower than the planet.
%                           ** Rotation is +delta.
%   * SUNLIT and FAST -- spacecraft flies by the planet on the sunlit side,
%                         and the spacecraft is faster than the planet.
%                           ** Rotation is -delta.
%   * DARKSIDE and SLOW -- spacecraft flies by the planet on the dark side,
%                         and the spacecraft is slower than the planet.
%                           ** Rotation is -delta.
%   * DARKSIDE and FAST -- spacecraft flies by the planet on the dark side,
%                         and the spacecraft is faster than the planet.
%                           ** Rotation is +delta.
%                         
switch lower(side)
  case 'sunlit'
    % SUNLIT and SLOW
    if( vB>vTB )
      phi2 = phi1 + delta;
    else
      % SUNLIT and FAST
      phi2 = phi1 - delta;
    end
  case {'darkside','dark'}
    % DARKSIDE and SLOW
    if( vB>vTB )
      phi2 = phi1 - delta;
    else
      % DARKSIDE and FAST
      phi2 = phi1 + delta;
    end
end
      
% now that we know the exit angle "phi2", we can compute the relative
% velocity vector when it exits the SOI:
vInfVector2 = vinf*[cos(phi2); sin(phi2)];

%% Compute the NEW spacecraft velocity vector and orbital elements "a", "e"
velVectorOfSpacecraftNewVS = velVectorOfPlanet + vInfVector2;

% the transverse component is in the positive "V" direction
vTrans = velVectorOfSpacecraftNewVS(1);
% the radial component is in the negative "S" direction
vRad = -velVectorOfSpacecraftNewVS(2);

% The delta-v:
dV = velVectorOfSpacecraftNewVS - velVectorOfSpacecraft;

% Given radial and transverse components of velocity, and "r", compute: a, e
h = rB*vTrans;
eSinTh = vRad*h/muS;
eCosTh = vTrans*h/muS - 1;
trueAnom = atan2( eSinTh, eCosTh );
ecc = eSinTh/sin(trueAnom);
sma = muS/(2*muS/rB-vTrans^2-vRad^2);

%% Plot the velocity vectors...
figure('name','Velocity Vectors')
quiver(0,0,velVectorOfPlanet(1),velVectorOfPlanet(2),1,'b')
hold on, grid on, axis equal
quiver(0,0,velVectorOfSpacecraft(1),velVectorOfSpacecraft(2),1,'g')
quiver(0,0,velVectorOfSpacecraftNewVS(1),velVectorOfSpacecraftNewVS(2),1,'r')
quiver(velVectorOfPlanet(1),velVectorOfPlanet(2),vInfVector1(1),vInfVector1(2),1,'c')
quiver(velVectorOfPlanet(1),velVectorOfPlanet(2),vInfVector2(1),vInfVector2(2),1,'m')
legend('Helio. Vel. of Panet B','Helio. Vel. of S/C Before Flyby',...
  'Helio. Vel. of S/C After Flyby','Rel. Vel. of S/C Before Flyby',...
  'Rel. Vel. of S/C After Flyby')

%% We could use the "Flyby" function to compute and plot stuff too...
rPH = [rB;0];
vPH = [0;vB];
vSH0 = [0;vTB];
d = Flyby( muS, muP, rPH, vPH, vSH0, rp, side, radP, sprintf('%s Flyby',planetB), 1 );

%% Plot all the orbits... 
%   * planet A's orbit
%   * planet B's orbit
%   * the Hohmann transfer ellipse
%   * the NEW orbit after flyby at planet B

figure('name','Orbits')
thx = linspace(0,2*pi,300);
plot(rA*cos(thx),rA*sin(thx),'b'), hold on, grid on, zoom on, axis equal
plot(rB*cos(thx),rB*sin(thx),'r')

% Hohmann transfer orbit
eT = abs(rB-rA)/(rB+rA);
aT = .5*(rA+rB);
[rT,xT,yT]=PerifocalOrbit( aT, eT, thx );
plot(xT,yT,'k')

% New orbit after flyby
r0 = [-rB;eps;eps];
v0 = [ [0 1;-1 0]*velVectorOfSpacecraftNewVS; 0]
[a,i,W,w,e,th0] = OrbitalElementsFromRV(r0,v0,muS)
posN = RVFromCOE(a,i,W,w,e,thx,muS);
plot(posN(1,:),posN(2,:),'g')
legend(planetA,planetB,'Hohmann Transfer','New Orbit After Flyby')

