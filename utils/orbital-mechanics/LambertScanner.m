function [dV1,dV2,tEnd,v1r,elT] = LambertScanner( r1,v1, r2,v2, startTimes, tOFs, mu )

% Call LambertSolver in a loop with an array of start times and an array of
% times of flight. For each one, compute the initial and final delta-v.
% The final delta-v is only needed if you want to rendezvous at the second point.
%
% The goal here is to find combinations of start-time and TOF that lead to
% smaller delta-v. 
%
% Note: This algorithm does not account for fly-by maneuvers explicitly.
%       However, you can use this in a post-flyby analysis. Just set v1 to
%       be the final heliocentric velocity of the spacecraft after it does
%       the flyby. ** Important: Only use startTimes=0 in this case. **
%
%   Inputs:
%     r1          Position of starting orbit at t=0
%     v1          Velocity of starting orbit at t=0
%     r2          Position of target orbit at t=0
%     v2          Velocity of target orbit at t=0
%     startTimes  Array of start times      (sec)
%     tOFs        Array of times of flight  (sec)
%     mu          Gravitational constant of central body
% 
%   Outputs:
%     dV1         Delta-v matrix from first orbit to transfer  (km/s)
%     dV2         Delta-v matrix from transfer to target orbit (km/s)
%     tEnd        Matrix of end-times, when the target orbit is reached (sec)
%     v1r         Cell array of required initial velocities (for transfer orbit)

if( nargin < 1 )

  echo on
  muS = 132.7e9;

  ss = SolarSystem(0);
  marsEl = ss(4).el;
  jupEl  = ss(5).el;

  % force zero inclination
  marsEl(2) = 0;
  jupEl(2)  = 0;

  % position and velocity of two planets at the same epoch
  [r1,v1] = RVFromCOE(marsEl,muS);
  [r2,v2] = RVFromCOE(jupEl, muS);
  
  startTimes  = (  0 : 30 : 1800)*86400;
  tOFs        = (200 : 60 : 3600)*86400;

  [dV1,dV2,tEnd] = LambertScanner( r1,v1, r2,v2, startTimes, tOFs, muS );
  echo off
  
  
  figure, surf(startTimes/86400/365,tOFs/86400/365,dV1')
  title('Departure Delta-V')
  xlabel('Start Time (yrs)')
  ylabel('Time of Flight (yrs)')
  zlabel('Delta-V (km/s)')

  figure, surf(startTimes/86400/365,tOFs/86400/365,dV2')
  title('Rendezvous Delta-V')
  xlabel('Start Time (yrs)')
  ylabel('Time of Flight (yrs)')
  zlabel('Delta-V (km/s)')
  
  return
end
  
nST = length(startTimes);
nTF = length(tOFs);

dV1 = zeros(nST,nTF);
dV2 = zeros(nST,nTF);
tEnd = zeros(nST,nTF);
v1r = cell(nST,nTF);
elT = v1r;

% loop over start times
for i=1:nST
  
  % start time (elapsed from epoch where init pos/vel are defined)
  t0 = startTimes(i);
  
  % propagate planet 1 up to the start-time
  [r1s,v1s] = RVAtTFromR0V0( r1, v1, t0, mu );
    
  % loop over times of flight
  for j=1:nTF
    
    % time of flight to use for interplanetary transfer maneuver
    tOF = tOFs(j);
    
    tEnd(i,j) = t0+tOF;
    
    % propagate planet 2 into future by start-time + time of flight
    [r2f,v2f] = RVAtTFromR0V0( r2, v2, t0 + tOF, mu );
    
    % solve Lambert's problem
    % Given
    %   rP1Now    = position of planet 1 (and spacecraft) now       (at t = 0)
    %   rP2Later  = position of planet 2 later                      (at t = tOF)
    %   vIPT1     = velocity of interplanetary transfer at planet 1 (at t = 0)
    %   vIPT2     = velocity of interplanetary transfer at planet 2 (at t = tOF)
    [vIPT1,vIPT2] = LambertSolver( r1s, r2f, tOF, mu, 'prograde' );
    
    % compute delta-vs
    dV1(i,j) = norm(vIPT1-v1s);
    dV2(i,j) = norm(vIPT2-v2f);
    
    v1r{i,j} = vIPT1;
    [aa,ii,WW,ww,ee,th] = OrbitalElementsFromRV(r1s,vIPT1,mu);
    elT{i,j} = [aa,ii,WW,ww,ee,th];
    
  end
  
end


