%% How to plot an orbit
%
% Here I describe 3 different ways to plot an orbit, starting with the
% position and velocity vector of the orbit at a particular time.
%
% Method 1:
% 1. Start with position and velocity vectors at time t0:  [r0, v0]
%
% 2. Compute the orbital elements: [a, i, W, w, e, th0]
%     - th0 is the true anomlay at time t0
%     - the true anomaly of the orbit will change over time, of course
%     - the other 5 elements are all fixed for this orbit.
%     - You can use the provided "OrbitalElementsFromRV.m" or use your own
%       function.
%
% 3. Create an array of true anomalies that are valid for this orbit
%     - For an elliptical orbit (e<1), the range [0, 2*pi] is sufficient.
%     - For a hyperbolic orbit (e>1), you have to take care to define the
%       true anomaly over its valid range. First compute the bounds:
%         thetaMin = -acos(-1/e);
%         thetaMax = acos(-1/e);
%       Then, define an array of values between these bounds.  Be careful
%       not to INCLUDE the bounds, because r will go to inifinity at those
%       values of theta. For example:
%         theta = linspace(thetaMin+0.01, thetaMax-0.01);
%
% 4. Compute [x,y] coordinates for the orbit in the perifocal frame over
%    the range of true anomalies you defined.
%
% 5. Compute the rotation matrix "M" that rotates from the perifocal frame
%     to the geocentric/equatorial frame (also called Earth-Centered
%     Inertial, or "ECI" for short, for Earth-centered orbits).
%
% 6. Compute [xI,yI,zI] coordinates for the orbit at all true anomaly
%    values. Do this by rotating each [x,y,0] vector from the perifocal
%    frame to the ECI frame, using the matrix "M".
%
% 7. Plot each [xI,yI,zI] coordinate in a 3D plot using the "plot3"
% command.
%
%
% Method 2:
% 1. Start with position and velocity vectors at time t0:  [r0, v0]
%
% 2. Define a "deltaTheta" vector that goes from 0 to 2*pi.
% 
% 2. Compute the Lagrange coefficients "f" and "g" at each value of dTheta.
%    To do this, implement equations 2.158 (a) thru (d) from the text.
%
% 3. Use the Lagrange coefficients to compute the position vector
%     [xI,yI,zI] associated with each dTheta.
%
% 4. Plot each [xI,yI,zI] coordinate in a 3D plot using the "plot3"
% command.
%
% 
% Method 3:
% 1. Start with position and velocity vectors at time t0:  [r0, v0]
%
% 2. Define a "right-hand-side" dynamics function for orbit dynamics with a
%    central body.
%
% 3. Define a time vector for integration. 
%     - For an elliptical orbit, you can just define it over one orbit 
%       period.
%     - For a hyperblic orbit, you will have to integrate backward in time
%       to get one portion of the trajectory, and forward in time to get
%       the other portion.
%     - Note that you would have to compute the orbital elements first in
%       order to know what type of orbit you have, and what the orbit 
%       period is!
%
% 4. Use ode45 or another integrator to compute the solution of [xI,yI,zI]
%     over time.
%
% 5. Plot each [xI,yI,zI] coordinate in a 3D plot using the "plot3"
% command.
%

%% Constants
mu = 398600.44;
Re = 6378.14;

%% Define the initial position and velocity
% Here are some examples. You can add more if you like!
% Pick one... comment out the others

% % a hyperbolic orbit
% r0 = [8778.7587135203;4779.78597603374;295.027919191783];
% v0 = [-5.23227615507787;9.54544661708752;1.04294540477248];
% 
% % a cicular orbit
% r0 = [8822.17134217377;4565.62475532981;1150.80988996769];
% v0 = [-2.93225354540026;5.07397112973843;2.34877630268796];

% an elliptical orbit
r0 = [-6299.71610945156;-2203.34734241343;2701.63606523828];
v0 = [3.23850017066731;-8.13490816482975;0.917075190747781];


%% Method 1: [r0,v0] -> orbital elements -> perifocal [x,y] -> ECI [x,y,z] across T.A.
[a,i,W,w,e,th0] = OrbitalElementsFromRV(r0,v0,mu);

if( e<1 )
  th = linspace(0,2*pi,300);
else
  thLim = acos(-1/e);
  th = linspace(-thLim+.01, thLim-.01, 300 );
end

% perifocal coordinates
r = a*(1-e^2)./(1+e*cos(th));
x = r.*cos(th);
y = r.*sin(th);

% ECI coordinates
Q = GeoEqToPerifocal(i,W,w);  % this rotates from ECI to perifocal
M = Q';              % the transpose rotates from perifocal to ECI

rECI = M*[x;y;zeros(size(x))];
xI = rECI(1,:);
yI = rECI(2,:);
zI = rECI(3,:);

% plot the orbit and add a sphere the size of earth
figure('name','METHOD 1'), plot3(xI,yI,zI,'-s')
grid on, axis equal, rotate3d on, xlabel('x'), ylabel('y'), zlabel('z')
hold on, [x,y,z]=sphere(50);
surf(x*Re,y*Re,z*Re,'facecolor','c','edgecolor','none')

%% Method 2: [r0,v0] -> Lagrange coeff's -> ECI [x,y,z] across T.A.
dTh = linspace(0,2*pi,300);
[f,g] = LagrangeCoeff(r0,v0,mu,dTh);
xI = f*r0(1) + g*v0(1);
yI = f*r0(2) + g*v0(2);
zI = f*r0(3) + g*v0(3);

% plot the orbit and add a sphere the size of earth
figure('name','METHOD 2'), plot3(xI,yI,zI,'-s')
grid on, axis equal, rotate3d on, xlabel('x'), ylabel('y'), zlabel('z')
hold on, [x,y,z]=sphere(50);
surf(x*Re,y*Re,z*Re,'facecolor','c','edgecolor','none')

% NOTE: if the orbit is hyperbolic, our definition of "dTh" will include 
% some true anomaly values that are beyond the asymptotic limits.
% This will cause result in us plotting a "reflection" of the actual orbit. 

% Here is a clever way to find those regions...
% Check if "f" makes a VERY large jump. 
% If it is hyperbolic, it will tend towards pos. or neg. infinity along
% the real trajectory, as dTh gets bigger.  
% When it switches from -inf to +inf, for example... That is
% when we cross the true anomaly bound.
if( max(abs(diff(f))) > 100 )
  
  k = find(abs(diff(f))>0.8*max(abs(diff(f))));
  real1 = 1:k(1);
  real2 = k(2)+1:length(dTh);
  
  % plot the orbit
  figure('name','METHOD 2 -- CHECKING FOR HYPERBOLIC'),
  plot3(xI(real1),yI(real1),zI(real1),'b-s'), hold on
  plot3(xI(real2),yI(real2),zI(real2),'g-s')
  
  % for illustrative purposes, plot the reflection
  ref = k(1)+1 : k(2);
  plot3(xI(ref),yI(ref),zI(ref),'r-o')
  legend('Real part future','Real part past','Reflection')
  
  % add a sphere the size of earth
  grid on, axis equal, rotate3d on, xlabel('x'), ylabel('y'), zlabel('z')
  hold on, [x,y,z]=sphere(50);
  surf(x*Re,y*Re,z*Re,'facecolor','c','edgecolor','none')
  
else
  
  
end

%% Method 3: [r0,v0] -> integrate dynamics -> ECI [x,y,z] across time

% "right-hand-side" dynamics function
% We will create "rhs" below. It is a function handle. 
%   - it accepts time "t" and state vector "x" as inputs.
%   - "rhs" returns d/dt(x), the time-rate-of-change of the state vector.
%   - The state vector is: [r; v]
%   - Therefore, the state derivative is: [d/dt(r); d/dt(v)]
%   - This is equivalent to:              [    v;      a   ]
%   - where "v" is the velocity, which is x(4:6)
%   - and "a" is the gravitational acceleration... -mu/r^2 in the "r" direction.
%   - Note that the time "t" is actually not needed, but we use it anyway 
%     because all of the ordinary differential equation (ODE) integrators
%     in Matlab expect a function handle that accepts (t,x) as inputs.
rhs = @(t,x) [x(4:6); ...
  -mu*x(1:3)/norm(x(1:3))^3];

% define a time vector...
% it would help to know the orbital elements so we can pick a suitable time
% vector...
[a,i,W,w,e,th0] = OrbitalElementsFromRV(r0,v0,mu);
if( e<1 )
  T = 2*pi/sqrt(mu/a^3);
  time = linspace(0,T,300);
  
  % integrate using ode45
  [t,xout] = ode45(rhs,time,[r0;v0],odeset('abstol',1e-12,'reltol',1e-12));
  xI = xout(:,1);
  yI = xout(:,2);
  zI = xout(:,3);

else
  % its hard to pick a good time frame with a hyperbolic orbit...
  % so just go back 12 hours and forward 12 hours
  time1 = linspace(0,-12,300)*86400;
  time2 = linspace(0,12,300)*86400;
  
  % integrate using ode45 backward in time
  [~,xout] = ode45(rhs,time1,[r0;v0],odeset('abstol',1e-12,'reltol',1e-12));
  xI1 = fliplr(xout(:,1)'); % "fliplr" flips the vector left/right, so it goes from -12 hrs to t0
  yI1 = fliplr(xout(:,2)');
  zI1 = fliplr(xout(:,3)');
  
  % integrate using ode45 forward in time
  [~,xout] = ode45(rhs,time2,[r0;v0],odeset('abstol',1e-12,'reltol',1e-12));
  xI2 = xout(:,1)';
  yI2 = xout(:,2)';
  zI2 = xout(:,3)';
  
  % stick the two segments together
  xI = [xI1,xI2];
  yI = [yI1,yI2];
  zI = [zI1,zI2];
  
end

% plot the orbit and add a sphere the size of earth
figure('name','METHOD 3'), plot3(xI,yI,zI,'-s')
grid on, axis equal, rotate3d on, xlabel('x'), ylabel('y'), zlabel('z')
hold on, [x,y,z]=sphere(50);
surf(x*Re,y*Re,z*Re,'facecolor','c','edgecolor','none')



