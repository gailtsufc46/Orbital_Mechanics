function [v1,v2,dTheta] = LambertSolver( r1, r2, TOF, mu, dir )

% Solve Lambert's problem.
%
%  Given two position vectors, r1 and r2, and the time of flight between
%  them, TOF, find the corresponding Keplerian orbit.
%

%% Demo
if( nargin<1 )
  mu = 398600.44;
  
  coe1 = [8000,pi/4,0,pi/3,.01,pi/6];
  coe2 = [12000,-pi/6,pi/4,2*pi/3,.2,pi/2];
  r1 = RVFromCOE(coe1,mu);
  r2 = RVFromCOE(coe2,mu);
  TOF = 12000;
  dir = 'pro';
  [v1,v2,dTheta] = LambertSolver( r1, r2, TOF, mu, dir );
  
  [fig,h1] = Plot3DOrbit(coe1,mu,0); 
  [~,h2] = Plot3DOrbit(coe2,mu,0,[],fig);
  
  coe3 = OrbitalElementsFromRV(r1,v1,mu);
  [~,h3] = Plot3DOrbit(coe3,mu,0,TOF,fig);
  disp(coe3)
  
  set(h1.orbit_plane,'facecolor','g','facealpha',.1)
  set(h2.orbit_plane,'facecolor','r','facealpha',.1)
  set(h3.orbit_plane,'facecolor','y','facealpha',.1)
  set(h1.orbit,'linewidth',2,'color',[0 .5 0])
  set(h2.orbit,'linewidth',2,'color',[1 0 0])
  set(h3.orbit,'linewidth',2,'color','k','linestyle','--')
  set(h3.arc,'color','m');
  set(h3.arc_start,'color','m','marker','.')
  set(h3.arc_stop,'color','m','marker','.')
  set([h1.x,h1.y,h1.z,h2.x,h2.y,h2.z,h3.x,h3.y,h3.z],'linewidth',1)
  
  return
end
  
%% Input checking

% direction
if( nargin<5 )
  dir = 'pro';
end

% Check validity of "dir" input
if( ~ischar(dir) )
  error('Input for direction "dir" must be a string, either ''pro'' or ''retro''.');
end

% mu
if( nargin<4 )
  mu = 398600.44;
  warning('Using Earth gravitational constant, mu = %f.',mu);
end

% TOF
if( TOF<=0 )
  error('The time of flight TOF must be >0.')
end

%% Calculate the position vector magnitudes
r1m = sqrt(r1'*r1);
r2m = sqrt(r2'*r2);

%% Calculate delta-theta using Equation 5.26
r1CrossR2 = cross(r1,r2);
arg = r1'*r2/r1m/r2m;
if( arg>1 )
  angle = 2*pi;
elseif( arg<-1 )
  angle = pi;
else
  angle = acos( arg );
end
switch lower(dir)
  case {'pro','prograde'}
    if( r1CrossR2(3) >= 0 )
      dTheta = angle;
    else
      dTheta = 2*pi-angle;
    end
  case {'retro','retrograde'}
    if( r1CrossR2(3) < 0 )
      dTheta = angle;
    else
      dTheta = 2*pi-angle;
    end
  otherwise
    error('Unrecognized direction. Use either ''pro'' or ''retro''.')
end

if( abs( abs(dTheta)-2*pi ) < 1e-8 )
  % this is a rectilinear orbit... not supported.
  v1 = nan;
  v2 = nan;
  return;
end

%% Calculate "A" in Equation 5.35
A = sin(dTheta)*sqrt(r1m*r2m/(1-cos(dTheta)));

% Inline function handles, all functions of "z"
yf  = @(z) r1m+r2m+A*( z.*stumpS(z)-1 )./sqrt(stumpC(z));
Ff  = @(z) ( yf(z)./stumpC(z) ).^1.5 .* stumpS(z) + A*sqrt(yf(z)) - sqrt(mu)*TOF;
dFf = @(z) LambertFPrimeOfZ( r1m, r2m, dTheta, z );

%% Solve for z ...

% Initial guess for z...
z0s = FindInitialZGuessForLambert(yf,Ff);
  
% consider each approx. crossing and solve for exact z value at each...
% terminate as soon as we have a good answer.
success = 0;
for i=1:length(z0s)
  [zs,success] = NewtonRhapsonSolver( Ff, dFf, z0s(i), 1e-6 );
  if( success )
    break;
  end
end
if( ~success )
  warning('Newton method did not find a solution in LambertSolver.')
  v1=nan;
  v2=nan;
  return;
end


%% Calculate "y" in Equation 5.38
%y = yf(zs);

%% Calculate the Lagrange Coefficients using Equation 5.46
[f,g,~,gdot] = LagrangeCoeffZ(r1m,r2m,dTheta,zs,mu);

%% Calculate v1 and v2 from Equation 5.28 and 5.29
v1 = 1/g*(r2-f*r1);
v2 = 1/g*(gdot*r2-r1);

%% Check result...

% Calculate the COE from [r1,v1]
[a1,i1,W1,w1,e1,th1] = OrbitalElementsFromRV(r1,v1,mu);

% Calculate the COE from [r2,v2]
[a2,i2,W2,w2,e2,th2] = OrbitalElementsFromRV(r2,v2,mu);

if( norm([1,cos([i1,W1,w1]),e1]-[a2/a1,cos([i2,W2,w2]),e2])>1e-8 )
  warning('Orbital elements from [r1,v1] do not match with those from [r2,v2].')
  disp([a1,i1,W1,w1,e1,th1; a2,i2,W2,w2,e2,th2]')
end

function z0 = FindInitialZGuessForLambertxxxxx( yf, Ff )

% A second attempt at a more robust way of finding a GOOD initial guess for
% the universal variable "z" where F(z)=0

zmin = -1e3;
zmax = 1e3;

% first find any y(z)=0 crossings, if they exists
zz = linspace(zmin,zmax,1e4);
yz = yf(zz);
y1 = yz(1:end-1);
y2 = yz(2:end);
kyz = find(sign(y1).*sign(y2)<0);

if( isempty(kyz) )
  %disp('No y(z)=0 crossing found...')
  zz = linspace(zmin,zmax,1e4);
  z0 = FindApproxFZero( Ff, zz );
  
else
  %disp(sprintf('%d y(z)=0 crossings found...',length(kyz)))
  z0 = [];
  for i=1:length(kyz)
    
    zy0 = fzero(yf,kyz(i));
    zz = linspace(zy0,zy0+zmax,1e4);
    z0i = FindApproxFZero( Ff, zz );
    if( ~isempty(z0i) )
      z0 = [z0, z0i];
    end
    
  end
end

%disp(z0)

function z0 = FindApproxFZero( Ff, zz )

% look over a range of z values to find approx. where F(z)=0

f = Ff(zz);                  % F(z) where y(z) > 0

% find approx. F(z)=0 crossing
f1 = f(1:end-1);
f2 = f(2:end);
ks = find(sign(f1).*sign(f2)<0);
z0 = zz(ks);

function Sz = SolveStumpffS(z)
%Used to determine stumpff equations S(z)
%for use in orbital equations (chp 3) and Lambert problem (chp 5)
%Input:
%   z - initial guess (updates)
%       *note: z = alpha*chi, alpha = 1/a, a = semimajor axis, 
%               chi = universal variable
%Outputs:
%   S(z) = Stemppf function

if z>0
    Sz = ( sqrt(z) - sin(sqrt(z)) )/(z^(3/2));
elseif z<0
    Sz = ( sinh(sqrt(-z)) - sqrt(-z) )/((-z)^(3/2));
else %z = 0
    Sz = 1/6;
end




