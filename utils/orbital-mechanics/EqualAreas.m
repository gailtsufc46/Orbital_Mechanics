function EqualAreas( a, e, th1, th2, dT, mu )
%
% Plot equal areas over equal times
%
% Inputs
%   a       Semi major axis (km)
%   e       Eccentricity (0 <= e < 1)
%   th1     True anomaly value where first area sweep begins (rad)
%   th2     True anomaly value where second area sweep begins (rad)
%   dT      The time period to sweep (sec)
%   mu      Gravitational parameter (km^3/s^2)

% demo (if no inputs are given)
if( nargin<1 )
  a = 8000;
  e = 0.4;
  th1 = 0;
  th2 = pi;
  dT = 600;
  mu = 398600.44;
end

if( nargin<6 )
  mu = 398600.44;
end

% true anomaly
th = linspace(0,2*pi,1e3);

% position
r = a*(1-e^2)./(1+e*cos(th));
x = r.*cos(th);
y = r.*sin(th);

% plot orbit
figure,
subplot(211)
plot(x,y), hold on, axis equal
subplot(212)
plot(x,y), hold on, axis equal

% start time for first area sweep
t1 = TimeSincePeriapsis(th1,a,e,mu);
% end time for first area
t1f = t1+dT;
% true anomaly at end of first area sweep
th1f = TrueAnomFromTime( t1f,a,e,mu );
if( th1f<th1 )
  th1f = th1f+2*pi;
end
th1x = linspace(th1,th1f,1e4);
r1x  = a*(1-e^2)./(1+e*cos(th1x));
subplot(211)
area1 = trapz(th1x,0.5*r1x.^2);
fill([0,r1x.*cos(th1x)],[0,r1x.*sin(th1x)],'r')
fprintf(1,'Numerical calculation of area 1: %f\n',area1)

% start time for second area sweep
t2 = TimeSincePeriapsis(th2,a,e,mu);
% end time for second area
t2f = t2+dT;
% true anomaly at end of second area sweep
th2f = TrueAnomFromTime( t2f,a,e,mu );
if( th2f<th2 )
  th2f = th2f+2*pi;
end
th2x = linspace(th2,th2f,1e4);
r2x  = a*(1-e^2)./(1+e*cos(th2x));
subplot(212)
area2 = trapz(th2x,0.5*r2x.^2);
fill([0,r2x.*cos(th2x)],[0,r2x.*sin(th2x)],'c')
fprintf(1,'Numerical calculation of area 2: %f\n',area2)

% analytically: dA/dt = h/2 (eq. 2.32)
h = sqrt(a*mu*(1-e^2));
dA = dT*h/2;
fprintf(1,'Analytical result for the area:  %f\n',dA)


