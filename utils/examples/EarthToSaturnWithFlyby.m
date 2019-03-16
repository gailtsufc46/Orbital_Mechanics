%% Example of planning a sequence of interplanetary maneuvers.
%
% In this case, we want to reach Saturn using as little Delta-V as
% possible.
%
% We will consider two possible trajectories:
%
%   1. a direct transfer from Earth to Saturn
%   2. a transfer to Jupiter, flyby at Jupiter, then rendezvous at Saturn
% 
% Assumptions:
%
% * We assume Keplerian orbits for all planets.
% * All planets are assumed to lie in the ecliptic plane (zero inclination).
% * We assume the earliest possible start date is May 1, 2018.
% * We will neglect the transfer in/out of parking orbits.

%% Constants 
mu = 1.327e11;              % grav. parameter of the sun
year = 2018;
month = 5;
day = 1;
jD0 = J0(year,month,day);   % Julian date for our earliest start date

% Use "SolarSystem" to get the planet radius and planet mu values 
s = SolarSystem([2018,5,1],0); 

muE = 398600.44;
muJup = s(5).earthMass*muE;
muSat = s(6).earthMass*muE;

radE   = s(3).radius;
radJup = s(5).radius;
radSat = s(6).radius;

soiJup = s(5).rSOI;

%% Julian day to calendar date function...
JDay2Date = @(x) datevec(x-1721058.5);

%% Method 1 - Direct Earth-to-Saturn Transfer
%
% We will try to find a good start date by first planning a Hohmann transfer. 
% This approach will assume circular orbits for Earth and Saturn. 
% 
% We will then consider different start dates and different durations, using
% the ACTUAL eccentric orbits for both planets.
% 
% So... we will have 2 independent variables:
% 
%     * jD1: Start Date
%     * tT:  Transfer time (duration)


% First get the radius and initial phase angle for both planets...
[pos1,~,~,coe1] = PlanetData(3, year, month, day, 0, 0, 0, mu);
[pos2,~,~,coe2] = PlanetData(6, year, month, day, 0, 0, 0, mu);
% use the semi-major axis for each planet's circular orbit radius
r1 = coe1(1); 
r2 = coe2(1);
% look at the position in the ecliptic plane to compute the initial phase
% angle
th10 = atan2(pos1(2),pos1(1));
th20 = atan2(pos2(2),pos2(1));
% Compute the delta-vs
hData = OrbMvrHohmann(r1,r2,r1,r2,mu);
% Compute the wait time, transfer time, and synodic period
[tW,tT,tSyn] = HohmannTransferAnimation( r1, r2, th10, th20, mu );

% If we had circular orbits, we would want to do a Hohmann transfer
% starting at day: jD0 + tW/86400
jDStartH = jD0 + tW/86400;

% But we have slightly eccentric orbits, so let's consider a range of start
% dates NEAR this date, and a range of transfer times NEAR the Hohmann
% transfer time...
jD1Vec = jDStartH + (-20 : 2 : 20); % e.g. from 50 days before to 50 days after
tTVec  = tT + (-200 : 2 : 300)*86400; % e.g. from 50 days shorter to 50 days longer

% make big empty matrices for first/second delta-v's
nST = length(jD1Vec); % number of different start times
nTT = length(tTVec); % number of different transfer times
dV1 = zeros(nST,nTT);
dV2 = zeros(nST,nTT);

% For each combination of start date and transfer time ... use Lambert to find
% the transfer orbit.
for i=1:length(jD1Vec)
  
  % start date
  date1 = JDay2Date(jD1Vec(i));
  y1 = date1(1);
  m1 = date1(2);
  d1 = date1(3);
  
  for j=1:length(tTVec)
    % position of Earth at start date
    [r1,v1] = PlanetData(3,y1,m1,d1,0,0,0,mu);
    
    % TOF, end date
    TOF = tTVec(j);
    date2 = JDay2Date( jD1Vec(i) + tTVec(j)/86400 );
    y2 = date2(1);
    m2 = date2(2);
    d2 = date2(3);

    % position of Saturn at start date
    [r2,v2] = PlanetData(6,y2,m2,d2,0,0,0,mu);
    
    % compute transfer maneuver
    [v1T,v2T] = LambertSolver(r1,r2,TOF,mu,'pro');
    
    % compute and record delta-v's
    dV1(i,j) = norm(v1T-v1);
    dV2(i,j) = norm(v2T-v2);
  end
  
end
dVTot = dV1+dV2;

% Plot the results of the Direct Earth-Saturn analysis...
figure, plot(jD1Vec-jDStartH,dVTot), grid on, xlabel('Start day variation'), ylabel('DV')
figure, surf(jD1Vec-jDStartH,(tTVec-tT)/86400,dVTot','edgecolor','none'), grid on, rotate3d on
xlabel('Start day variation (days)'), ylabel('Transfer time variation (days)'), zlabel('DV')
hold on, plot3(0,0,hData.dVTot,'r.','markersize',30)
shading interp
[v,row]=min(dVTot); [v,col]=min(v); row=row(col);
plot3(jD1Vec(row)-jDStartH,(tTVec(col)-tT)/86400,dVTot(row,col),'g.','markersize',30)
legend('Lambert Solver DVs','Hohmann Solution','Min DV from Lambert','location','northeast')

dVTot_method1 = min(dVTot(:));
days_method1  = tTVec(col)/86400;

jD1_1 = jD1Vec(row);
jD2_1 = jD1_1 + days_method1;


%% Method 2 -  Earth Transfer to Jupiter, flyby at Jupiter, then rendezvous at Saturn
% 
% Step 1 - Transfer to Jupiter...
% +++++++++++++++++++++++++++++++
% (similar to above Earth-Saturn analysis)
% 
% We will try to find a good start date by first planning a Hohmann transfer
% to Jupiter.
% This approach will assume circular orbits for Earth and Jupiter. 
% 
% We will then consider different start dates and different durations, using
% the ACTUAL eccentric orbits for both planets.
% 
% So... we will have 2 independent variables:
%     * jD1:  Start Date
%     * tT1:  Transfer time (duration) from Earth to Jupiter
% 
% Using Lambert, we will compute v1T and v2T. Here, "v2T" is the velocity of
% our spacecraft in its transfer orbit when we arrive at Jupiter. This is in
% the heliocentric frame. If we subtract Jupiter's velocity at this time
% (v2J), we will get the incoming hyperbolic excess velocity vector,
% "vInfVec1".
% 
% From this analysis, we will get our first required delta-v, "dV1". This is
% applied at date "jD1" to leave Earth and head to Jupiter.
%     dV1 = norm( v1E - v1T );
% 
% Step 2 - Compute the next Transfer to Saturn
% ++++++++++++++++++++++++++++++++++++++++++++
% 
% We are at Jupiter (r2) at date jD2. We want to arrive at Saturn some time later.
% We get to choose the transfer time, tT2. 
% We will define the final position "r3" to be the position of Saturn at date
% jD3, where jD3 = jD2 + tT2/86400.
% We will use Lambert to compute the transfer orbit from r2 to r3, with a
% time of flight of tT2.
% 
% So, for this step, we have just one variable:
%     * tT2:  Transfer time (duration) from Jupiter to Saturn
% 
% Using Lambert, we will compute v2Tp and v3T, where "v2p" is the required
% new velocity (in the heliocentric frame) of our spacecraft post-flyby.
% 
% This analysis also tells us "dV3", the final delta-v to be applied at
% Saturn to rendezvous with that planet. It is equal to:
%     dV3 = norm( v3S - v3T );
% 
% Step 3 - Compute the best possible flyby at Jupiter
% +++++++++++++++++++++++++++++++++++++++++++++++++++
% Given "vInfVec1" that we found from Step 1, and "v2req" that we found from
% step 2... our goal in this step is to choose the best possible flyby
% maneuver. How do we define "best" here? Remember that the result of a flyby
% is a new heliocentric velocity. We would like our flyby maneuver to give us
% a new heliocentric velocity that is AS CLOSE AS POSSIBLE to the required
% heliocentric velocity, "v2req". 
% 
% For a flyby maneuver, we can choose the periapse distance of the hyperbolic
% orbit, "rp". We are constrained that rp > planetRadius. 
% 
% So, for this step, we have just one variable: 
%   
%   * rp:     Periapse distance of the hyperbolic orbit.
% 
% Using our flyby analysis, we compute the new heliocentric velocity of our
% spacecraft AFTER the flyby, and compare that to the required velocity to
% reach Saturn.
%   v2PostFlyby = v2J + vInfVec2;
%   dV2 = norm( v2PostFlyby - v2p );
% 
% Choose rp to give us the smallest dV2 possible.
% 
%  SUMMARY
% +++++++++++++++++++
% 
% * Do Step 1 just once to get a good Earth-Jupiter transfer. Output --> dV1
% * Do Steps 2 and 3 together. 
%  -   Pick a Jupiter-to-Saturn transfer time, "tT2". Output --> dV2, dV3
%  -  Compute total delta-v:  dVTot = dV1 + dV2 + dV3
%  -  Repeat with different transfer times "tT2" until dVTot is minimum.


%% Method 2 - Step 1 

% First get the radius and initial phase angle for both planets...
[pos1,~,~,coe1] = PlanetData(3, year, month, day, 0, 0, 0, mu);
[pos2,~,~,coe2] = PlanetData(5, year, month, day, 0, 0, 0, mu);
% use the semi-major axis for each planet's circular orbit radius
r1 = coe1(1); 
r2 = coe2(1);
% look at the position in the ecliptic plane to compute the initial phase
% angle
th10 = atan2(pos1(2),pos1(1));
th20 = atan2(pos2(2),pos2(1));
% Compute the delta-vs
hData = OrbMvrHohmann(r1,r2,r1,r2,mu);
% Compute the wait time, transfer time, and synodic period
[tW,tT,tSyn] = HohmannTransferAnimation( r1, r2, th10, th20, mu );

% If we had circular orbits, we would want to do a Hohmann transfer
% starting at day: jD0 + tW/86400
jDStartH = jD0 + tW/86400;

% But we have slightly eccentric orbits, so let's consider a range of start
% dates NEAR this date, and a range of transfer times NEAR the Hohmann
% transfer time...
jD1Vec = jDStartH + (-20 : 2 : 20); % e.g. from 50 days before to 50 days after
tTVec  = tT + (-200 : 2 : 200)*86400; % e.g. from 50 days shorter to 50 days longer

% make big empty matrices for first/second delta-v's
nST = length(jD1Vec); % number of different start times
nTT = length(tTVec); % number of different transfer times
dV1 = zeros(nST,nTT);
dV2 = zeros(nST,nTT);

% For each combination of start date and transfer time ... use Lambert to find
% the transfer orbit.
for i=1:length(jD1Vec)
  
  % start date
  date1 = JDay2Date(jD1Vec(i));
  y1 = date1(1);
  m1 = date1(2);
  d1 = date1(3);
  
  for j=1:length(tTVec)
    % position of Earth at start date
    [r1,v1] = PlanetData(3,y1,m1,d1,0,0,0,mu);
    
    % TOF, end date
    TOF = tTVec(j);
    date2 = JDay2Date( jD1Vec(i) + tTVec(j)/86400 );
    y2 = date2(1);
    m2 = date2(2);
    d2 = date2(3);

    % position of Earth at start date
    [r2,v2] = PlanetData(5,y2,m2,d2,0,0,0,mu);
    
    % compute transfer maneuver
    [v1T,v2T] = LambertSolver(r1,r2,TOF,mu,'pro');
    
    % compute and record delta-v's
    dV1(i,j) = norm(v1T-v1);
    dV2(i,j) = norm(v2T-v2);
  end
  
end
dVTot = dV1+dV2;

% Plot the results of the Direct Earth-Jupiter analysis...
figure, plot(jD1Vec-jDStartH,dVTot), grid on, xlabel('Start day variation'), ylabel('DV')
figure, surf(jD1Vec-jDStartH,(tTVec-tT)/86400,dVTot','edgecolor','none'), grid on, rotate3d on
xlabel('Start day variation (days)'), ylabel('Transfer time variation (days)'), zlabel('DV')
hold on, plot3(0,0,hData.dVTot,'r.','markersize',30)
shading interp
[v,row]=min(dVTot); [v,col]=min(v); row=row(col);
plot3(jD1Vec(row)-jDStartH,(tTVec(col)-tT)/86400,dVTot(row,col),'g.','markersize',30)
legend('Lambert Solver DVs','Hohmann Solution','Min DV from Lambert','location','northeast')

% choose best start date and transfer time:
jD1 = jD1Vec(row);
tT1 = tTVec(col);

% Dates jD1 and jD2
date1 = JDay2Date(jD1);
y1 = date1(1); m1 = date1(2); d1 = date1(3);

jD2 = jD1+tT1/86400;
date2 = JDay2Date(jD2);
y2 = date2(1); m2 = date2(2); d2 = date2(3);

% Recompute this transfer
[r1E,v1E] = PlanetData(3,y1,m1,d1,0,0,0,mu);
[r2J,v2J] = PlanetData(5,y2,m2,d2,0,0,0,mu);
[v1T,v2T] = LambertSolver(r1E,r2J,tT1,mu,'pro');

% dV1
dV1 = norm( v1T - v1E );

%% Method 2 - Steps 2,3

% Define a range of transfer times from Jupiter to Saturn
%tT2x = (5:2:15) * 365*86400;
tT2x = (12: .2 :14) * 365*86400;

% Define a set of periapse distances
rpx = linspace( radJup, radJup + 0.25*soiJup, 300 );  

figure('name','Delta-V 2...')
dV2 = zeros(length(tT2x),length(rpx));
dV3 = zeros(length(tT2x),length(rpx));
dV23 = zeros(length(tT2x),length(rpx));

for i=1:length(tT2x)
  
  % Pick the transfer time from Jupiter to Saturn
  tT2 = tT2x(i);
  
  % Date jD3
  jD3 = jD2 + tT2/86400;
  date3 = JDay2Date(jD3);
  y3 = date3(1); m3 = date3(2); d3 = date3(3);
  
  % Saturn at jD3
  [r3S,v3S] = PlanetData(6,y3,m3,d3,0,0,0,mu);
  
  % Lambert solution to transfer from Jupiter to Saturn
  [v2Tp,v3T] = LambertSolver(r2J,r3S,tT2,mu,'pro');
  
  for k=1:length(rpx)
    
    % Now the flyby analysis...
    fs = Flyby(mu,muJup,r2J,v2J,v2T,rpx(k),'sunlit',radJup,'',0);
    fd = Flyby(mu,muJup,r2J,v2J,v2T,rpx(k),'dark',radJup,'',0);
    
    % post-flyby heliocentric felocity...
    %v2PFB = [fd.vFH;0]; % NOTE: we could use either the darkside or sunlit side!
    v2PFB = [fs.vFH;0]; % NOTE: we could use either the darkside or sunlit side!
    
    % compute dV2...
    dV2(i,k) = norm( v2Tp - v2PFB );
    
    % compute dV3...
    dV3(i,k) = norm( v3S - v3T );
    
    % and the total for later analysis / plotting
    dV23(i,k)=dV2(i,k)+dV3(i,k);
    
  end
  
  % Plot the variation in dV2 across the range of periapse distances
  subplot(311)
  plot(rpx,dV2(i,:),'displayname',sprintf('%1.1f Years Transfer',tT2/86400/365))
  hold on

  subplot(312)
  plot(rpx,dV3(i,:),'displayname',sprintf('%1.1f Years Transfer',tT2/86400/365))
  hold on

  subplot(313)
  plot(rpx,dV23(i,:),'displayname',sprintf('%1.1f Years Transfer',tT2/86400/365))
  hold on

end
subplot(311), grid on
ylabel('DV2 (km/s)')
subplot(312), grid on
ylabel('DV3 (km/s)')
subplot(313), grid on
ylabel('DV2 + DV3 (km/s)')
xlabel('Periapse Distance, rp (km)')
legend

% Find the minimum
[v,rows]=min(dV23); [v,col]=min(v); row=rows(col);

% Date jD3
tT2opt = tT2x(row);
jD3 = jD2 + tT2opt/86400;
date3 = JDay2Date(jD3);
y3 = date3(1); m3 = date3(2); d3 = date3(3);

% Saturn at jD3
[r3S,v3S] = PlanetData(6,y3,m3,d3,0,0,0,mu);

% Corresponding transfer orbit
[v2Tp,v3T] = LambertSolver(r2J,r3S,tT2opt,mu,'pro'); 

%% So what?
%{
Compute the minimum total delta-v for method 2...

Is it smaller than the total transfer delta-v for method 1? If so, cool. 
How long does it take?

%}
dV2min = dV2(row,col);
dV3min = dV3(row,col);
dVTot_method2 = dV1 + dV2min + dV3min;
days_method2 = tT2opt/86400;

fprintf(1,'Method 1: Total Transfer Delta-V: %2.3f km/s\n',dVTot_method1);
fprintf(1,'Method 2: Total Transfer Delta-V: %2.3f km/s\n\n',dVTot_method2);
fprintf(1,'Method 1: Total Transfer Time: %4.1f days\n',days_method1);
fprintf(1,'Method 2: Total Transfer Time: %4.1f days\n',days_method2);

%% Let's record the 3 delta-v's experienced by the spacecraft

% The spcacecraft is initially at Earth, moving at Earth's velocity
[r0,v0] = PlanetData(3,year,month,day,0,0,0,mu);
sc(1).r0 = r0;
sc(1).v0 = v0;

% First delta-v is the first burn for the Earth-Jupiter Lambert Transfer,
% to depart Earth. (We ignore parking orbits)
sc(1).dV(:,1) = v1T-v1E; 
sc(1).dVJD(1) = jD1;

% Second TOTAL delta-v is the difference between "v2T" and "v2Tp", where:
%   v2T is the velocity of our Earth-Jupiter transfer orbit when arrive at Jupiter
%   v2Tp is the velocity of our Jupiter-Saturn transfer orbit when we depart Jupiter 
% So... this delta-v includes the DV from the flyby, plus whatever our
% propulsion system provides.
sc(1).dV(:,2) = v2Tp-v2T;
sc(1).dVJD(2) = jD2;

% Third delta-v is the final burn for the Jupiter-Saturn Lambert Transfer,
% to rendezvous with Saturn. (Again, we ignore parking orbits)
sc(1).dV(:,3) = v3S-v3T;
sc(1).dVJD(3) = jD3;

% Now animate the whole thing to see what's going on...
InterplanetaryAnimation([3 5 6],sc,jD0);

