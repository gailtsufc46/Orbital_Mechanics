%% Flyby to Intercept Example
%
% This example shows how to do a flyby maneuver followed by a delta-v (from
% your propulsion system) in order to intercept another planet.
% 
% Overall process:
%   + We approach planet 1 at time t = 0 with some initial heliocentric
%       velocity (vH1a).
%   + We do a flyby of planet 1. We can choose: 
%       - Dark side or Sunlit side
%       - Aiming radius (or equivalently, the periapsis distance "rP")
%   + At the end of the flyby, we have a new heliocentric velocity (vH1b).
%   + We use Lambert to solve for a DESIRED heliocentric velocity (vH1des). 
%       This is the velocity that we need to intercept with Planet 2 after a
%       specified time. Here, we can choose:
%       - Time of flight (tOF)
%   + The difference between our post-flyby velocity (vH1b) and our desired
%     velocity (vH1des) is the required delta-v that must be provided by our
%     propulsion system.
% 
%   + NOTE: Because SOI << R (to sun), we assume the spacecraft is located at 
%     the same point as the planet in the heliocentric frame.
% 
%   + NOTE: Because the time it takes to traverse the SOI is small compared
%     to the orbit period of the planets, we assume the planets do not move
%     during our flyby. This is equivalent to assuming that the delta-v 
%     resulting from our flyby is instantaneous. So its just like an
%     impulsive maneuver.
%

%% Constants
au = 149597871;
muS = 132.7e9;


%% Initial Conditions
ss = SolarSystem;

kp = 5;   % current planet index (3 is Earth, 4 is Mars, 5 is Jupiter)
kn = 6;   % planet index to visit next

el1 = ss(kp).el;
el2 = ss(kn).el;

% force zero inclination!!
el1(2) = 0;
el2(2) = 0;

% position and velocity of two planets at the same epoch
[rP1Now,vP1Now] = RVFromCOE(el1, muS);
[rP2Now,vP2Now] = RVFromCOE(el2, muS);


%% Examine the flyby occuring at Planet 1
muP   = ss(kp).mu;
radP  = ss(kp).radius;
rSOI  = ss(kp).rSOI;

% define an intial heliocentric velocity for the spacecraft ... 
% pick one randomly, but keep it fairly close to the planet's velocity
ra = rand*2*pi;
vH1a = vP1Now + 6*[cos(ra);sin(ra);0];

f_sunlit_close  = Flyby(muS,muP,rP1Now,vP1Now,vH1a, radP, 'sunlit', radP, 'Sunlit / Close');
f_sunlit_far    = Flyby(muS,muP,rP1Now,vP1Now,vH1a, rSOI, 'sunlit', radP, 'Sunlit / Far'  );
f_dark_close    = Flyby(muS,muP,rP1Now,vP1Now,vH1a, radP, 'dark',   radP, 'Dark Side / Close');
f_dark_far      = Flyby(muS,muP,rP1Now,vP1Now,vH1a, rSOI, 'dark',   radP, 'Dark Side / Far');

% two options for "vSCNow" - velocity of spacecraft now (right after flyby)
vSCNow_a = [f_sunlit_close.vFH;0];
vSCNow_b = [f_dark_close.vFH;0];


%% Use Lambert to examine DV requirement to intercept NEXT planet...
startTimes = 0;
%tOFs = 0+(18 : 2 : 365*3)*86400;
tOFs = 0+(230 : 5 : 365*3.5)*86400;
[dV1a,~,~,v1ra,elTa] = LambertScanner( rP1Now, vSCNow_a, rP2Now,vP2Now, 0, tOFs, muS );
[dV1b,~,~,v1rb,elTb] = LambertScanner( rP1Now, vSCNow_b, rP2Now,vP2Now, 0, tOFs, muS );

figure, plot(tOFs/86400/365,[dV1a;dV1b])
grid on, xlabel('TOF (yrs)'), ylabel('Delta-V (km/s)')
legend('Sunlit / Close','Dark side / Close')

%% solve Lambert's problem for a specific TOF
% Given
%   rP1Now    = position of planet 1 (and spacecraft) now       (at t = 0)
%   rP2Later  = position of planet 2 later                      (at t = tOF)
%   vIPT1     = velocity of interplanetary transfer at planet 1 (at t = 0)
%   vIPT2     = velocity of interplanetary transfer at planet 2 (at t = tOF)
[dVMin,kmin] = min(dV1a);
tOF = tOFs(kmin);
rP2Later = RVAtTFromR0V0( rP2Now, vP2Now, tOF, muS );
[vIPT1,vIPT2] = LambertSolver( rP1Now, rP2Later, tOF, muS, 'prograde' );

%% Velocity Vector Analysis...

PlotFlybyVector(f_dark_close);
title('Sunlit / Close')
for i=1:length(v1ra), 
  ph(i) = plot(v1ra{i}(1),v1ra{i}(2),'g.','handlevisibility','off','markersize',20); 
end
plot(v1ra{1}(1),v1ra{1}(2),'c*','displayname','Required New Vel. / Earliest TOF');
plot(v1ra{end}(1),v1ra{end}(2),'ms','displayname','Required New Vel. / Latest TOF');

% PlotFlybyVector(f_dark_close);
% title('Dark side / Close')
% for i=1:length(v1rb), 
%   ph(i) = plot(v1rb{i}(1),v1rb{i}(2),'g.'); 
% end
% plot(v1rb{1}(1),v1rb{1}(2),'c*');
% plot(v1rb{end}(1),v1rb{end}(2),'ms');



