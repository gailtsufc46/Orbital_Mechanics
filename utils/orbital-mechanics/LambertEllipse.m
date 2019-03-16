function d = LambertEllipse( a, c, r1r2, EF )

% Compute the two ellipse solutions for Lambert's theorem

% Lambert's theorem states that the time of flight between two points on a
% conic section (not just an ellispe) can be determined by JUST the
% following three pieces of information:
%   1. The semi major axis "a"
%   2. The chord length "c", which is just the distance between the points
%   3. The sum of the distances to each point "r1+r2" (from central body)
%
%   Inputs: 
%     a     Semi major axis
%     c     Chord (distance between r1 and r2)
%     r1r2  Sum of r1 and r2
%     EF    

% demo
if( nargin < 1 )
  a = 1.25;
  c = 1;
  r1r2 = 2;
  EF = 1.63;
end

% number of points for all curves
nn = 250;

% ellipse of the primary focus
aF = .5*(r1r2);
eF = c/(2*aF);

% ellipse of vacant focus
aFS = 2*a - aF;
eFS = c/(2*aFS);


% new plot with points p1 and p2 on horizontal line, x axis
x10 = 0; y10 = 0;
x20 = c; y20 = 0;
% figure
% plot(x1,y1,'b.','markersize',40,'linewidth',3), hold on
% plot(x2,y2,'r.','markersize',40,'linewidth',3)
% plot([x1 x2],[y1 y2],'k--')

Ex = linspace(0,2*pi,nn);
xFx = aF*cos(Ex)+c/2;
yFx = aF*sqrt(1-eF^2)*sin(Ex);
xFSx = aFS*cos(Ex)+c/2;
yFSx = aFS*sqrt(1-eFS^2)*sin(Ex);

% plot(xFx,yFx,'k',xFSx,yFSx,'m')
% axis equal

d.x10 = x10;
d.y10 = y10;
d.x20 = x20;
d.y20 = y20;
d.aF = aF;
d.eF = eF;
d.aFS = aFS;
d.eFS = eFS;
d.xFx = xFx;
d.yFx = yFx;
d.xFSx = xFSx;
d.yFSx = yFSx;

%% pick an arbitrary main focus 
xF = aF*cos(EF)+c/2;
yF = aF*sqrt(1-eF^2)*sin(EF);
% plot(xF,yF,'ks')
% plot([0 xF],[0 yF],'k')

% distance from each point to the main focus
dF1 = sqrt((xF-x10)^2+(yF-y10)^2);
dF2 = sqrt((xF-x20)^2+(yF-y20)^2);

% remaining distance from each point to OTHER empty focus
dFS1 = 2*a-dF1;
dFS2 = 2*a-dF2;

% draw circles with these radii around each point
th = linspace(0,2*pi,nn);
xF1c = x10+dFS1*cos(th); yF1c = y10+dFS1*sin(th);
xF2c = x20+dFS2*cos(th); yF2c = y20+dFS2*sin(th);
% plot(xF1c,yF1c,'b',xF2c,yF2c,'r')

d.xF = xF;
d.yF = yF;
d.th = th;
d.xF1c = xF1c;
d.yF1c = yF1c;
d.xF2c = xF2c;
d.yF2c = yF2c;
d.r1  = dF1;
d.r2  = dF2;
d.rS1 = dFS1;
d.rS2 = dFS2;

%% compute the exact possible locations of the empty focus
% its where the circles intersect... note that this intersection is also on
% the ellipse locus of points that we said the focus had to be on!

% define the origin at P1
xFS  = c/2 + (dFS1^2-dFS2^2)/(2*c);
yFS1 = +sqrt( ( (dFS1+dFS2)^2-c^2)*(c^2-(dFS1-dFS2)^2) ) / (2*c);
yFS2 = -sqrt( ( (dFS1+dFS2)^2-c^2)*(c^2-(dFS1-dFS2)^2) ) / (2*c);

% plot(xFS*[1 1],[yFS1 yFS2],'k*')

d.xFS = xFS;
d.yFS1 = yFS1;
d.yFS2 = yFS2;

%% compute the eccentricity
% distance between foci
dfoci1 = sqrt( (xFS-xF)^2 + (yFS1-yF)^2 );
dfoci2 = sqrt( (xFS-xF)^2 + (yFS2-yF)^2 );

% 2*a*e = distance between foci
e1 = dfoci1/(2*a);
e2 = dfoci2/(2*a);

d.e1 = e1;
d.e2 = e2;

%% consider each possible empty focus ... and compute the ellipse for it

Ep = linspace(0,2*pi,nn);

% orientation of each ellipse
alpha1 = atan2(yFS1-yF,xFS-xF);
alpha2 = atan2(yFS2-yF,xFS-xF);

% locus of points for ellipse 1 in perifocal frame
x1p = a*cos(Ep)-a*e1;
y1p = a*sqrt(1-e1^2)*sin(Ep);

% locus of points for ellipse 2 in perifocal frame
x2p = a*cos(Ep)-a*e2;
y2p = a*sqrt(1-e2^2)*sin(Ep);

% now transform back into our original frame
% we made P1 the origin...

% first ellipse
x1 = xF + x1p*cos(pi-alpha1) + y1p*sin(pi-alpha1);
y1 = yF - x1p*sin(pi-alpha1) + y1p*cos(pi-alpha1);

% second ellipse
x2 = xF + x2p*cos(pi-alpha2) + y2p*sin(pi-alpha2);
y2 = yF - x2p*sin(pi-alpha2) + y2p*cos(pi-alpha2);

% plot the ellipses
% plot(x1,y1,'g--','linewidth',2)
% plot(x2,y2,'c:','linewidth',2)

d.x1 = x1;
d.y1 = y1;
d.x2 = x2;
d.y2 = y2;



