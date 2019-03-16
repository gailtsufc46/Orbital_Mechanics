%% AEM 4301 - Compare a 1-burn Apse Rotation with a 2-burn Chase
%
% For an apse line rotation, the two orbits intersect at two points. Each
% of those points is an oportunity for a single-burn transfer maneuver.
%
% As an alternative, we use Lambert's method to plan a 2-burn maneuver
% that takes us from one point on the first orbit, to another point on the
% second orbit. 

%% Constants
mu = 398600.44;
Re = 6378.14;

%% Grid Spacing
dThetax = 4;      % angular spacing for Lambert search [deg]
dTx     = 0.025;  % time spacing for Lambert search [fraction of orbit 1 period]

%% Initial and Final orbits
a1 = 15000;
e1 = 0.5;
w1 = 0;     % orbit 1 argument of perigee

a2 = 15000;
e2 = 0.3;
w2 = 43.5*pi/180;  % orbit 2 argument of perigee

% "eta" is the angle between the two apse lines.
% Note that later we will solve for the true latitude (of each orbit) 
% where the two orbits intersect. It will also be the case that "eta" is
% equal to theta2 - theta1 (both at point "I" and at point "J").
eta = w2-w1;

%% Function handles to compute position and velocity vectors in common frame
posVecAtTA = @(h,e,w,theta) ... 
  h^2/mu/(1+e*cos(theta))*RotMat(w,3)*[cos(theta);sin(theta);0];

velVecAtTA = @(h,e,w,theta) ... 
  mu/h*RotMat(w,3)*[-sin(theta); e+cos(theta); 0];


%% Plot the orbits
np = 200;
theta = linspace(0,2*pi,np);
[r1,x1,y1] = PerifocalOrbit( a1, e1, theta );
[r2,x2,y2] = PerifocalOrbit( a2, e2, theta );
mat1 = RotMat(w1,3);
mat2 = RotMat(w2,3);
pos1 = mat1*[x1;y1;zeros(1,np)];
pos2 = mat2*[x2;y2;zeros(1,np)];

f1=figure; plot(pos1(1,:),pos1(2,:)), hold on, plot(pos2(1,:),pos2(2,:))
grid on, axis equal, zoom on
fill(Re*cos(theta),Re*sin(theta),'b')

%% Apse line rotation

% first compute h given a,e (for each orbit)
h1 = sqrt(a1*mu*(1-e1^2));
h2 = sqrt(a2*mu*(1-e2^2));

% follow the solution procedure
A = e1*h2^2 - e2*h1^2*cos(eta);
B = -e2*h1^2*sin(eta);
C = h1^2-h2^2;
phi = atan2(B,A);
th1i = phi + acos(C/A*cos(phi));
th1j = phi - acos(C/A*cos(phi));

pos1i = posVecAtTA(h1,e1,w1,th1i);
pos1j = posVecAtTA(h1,e1,w1,th1j);

plot(pos1i(1),pos1i(2),'ko','markersize',16,'linewidth',3)
plot(pos1j(1),pos1j(2),'ks','markersize',16,'linewidth',3)

% true anomaly values for orbit 2
th2i = th1i-eta;
th2j = th1j-eta;

% delta-v at point I
vel1i = velVecAtTA(h1,e1,w1,th1i);
vel2i = velVecAtTA(h2,e2,w2,th2i);
dVi   = norm(vel2i-vel1i);

% delta-v at point J
vel1j = velVecAtTA(h1,e1,w1,th1j);
vel2j = velVecAtTA(h2,e2,w2,th2j);
dVj   = norm(vel2j-vel1j);

fprintf(1,'\nAPSE LINE ROTATION:\n\tDelta-v at I: %f km/s\n\tDelta-v at J: %f km/s\n\n',...
  dVi,dVj);


%% Two Burn maneuver using Lambert Solver

% array of angles 
thx1 = (0:dThetax:(360-dThetax))*pi/180;
thx2 = (0:dThetax:(360-dThetax))*pi/180 + randn*pi/180;
npts1 = length(thx1);
npts2 = length(thx2);

% orbit period of first orbit
T1 = OrbPeriod(a1,mu);

% transfer times
ttimes = (.1 : dTx : 1.0)*T1;
nt = length(ttimes);

dVTot = ones(npts1,npts2)*inf;
bestTime = zeros(npts1,npts2);

fprintf(1,'Running %d Lambert solutions...\n',npts2*npts1*nt);

ts=tic;
deltaPerc = 5;
nextPercMarker = deltaPerc;
for i=1:npts1
  
  r1x = posVecAtTA( h1, e1, w1, thx1(i) );
  v1x = velVecAtTA( h1, e1, w1, thx1(i) );
  
  for j=1:npts2
   
    r2x = posVecAtTA( h2, e2, w2, thx2(j) );
    v2x = velVecAtTA( h2, e2, w2, thx2(j) );
    
    for k=1:nt
      
      [v1,v2] = LambertSolver(r1x,r2x,ttimes(k),mu,'pro');
      if( isnan(v1) )
        dVTmp = inf;
      else
        dVTmp = norm(v1-v1x) + norm(v2-v2x);
      end
      if( dVTmp < dVTot(i,j) )
        dVTot(i,j) = dVTmp;
        bestTime(i,j) = ttimes(k);
      end
      
    end
    
    percDone = i*j/npts1/npts2*100;
    if( percDone >= nextPercMarker )
      nextPercMarker = nextPercMarker + deltaPerc;
      if( percDone < 100 )
        remPerc = 100-percDone;
        estTimeLeft = toc/deltaPerc*remPerc;
        fprintf(1,'%1.1f %% complete. Approx. %2.1f seconds until finished.\n',percDone,estTimeLeft);
      end
      tic
    end
  end
end

% get rid of the "inf" values for plotting
actualMaxDV = max(dVTot(dVTot<inf));
dVTot(isinf(dVTot)) = actualMaxDV;

%% plot the delta-vs for all starting/ending points
f2 = figure; ss=surf(thx1*180/pi,thx2*180/pi,(dVTot'));
set(ss,'edgecolor','none'), shading interp, lighting phong
colormap colorcube
view(0,90)
xlabel('Orbit 1 True Anomaly (deg)')
ylabel('Orbit 2 True Anomaly (deg)')
zlabel('Total Delta-V (km/s)')
colorbar

%% compute and plot minimum delta-v transfer orbit
[tmp,rows] = min(dVTot);
[minDVTot,col] = min(tmp);
row = rows(col);
r1x = posVecAtTA( h1, e1, w1, thx1(row) );
r2x = posVecAtTA( h2, e2, w2, thx2(col) );
dT  = bestTime(row,col);
[v1,v2] = LambertSolver(r1x,r2x,dT,mu,'pro');

[rt,vt] = RVAtTFromR0V0(r1x,v1,linspace(0,dT,200),mu);

figure(f1)
hold on
plot(rt(1,:),rt(2,:),'g','linewidth',2)
plot(r1x(1),r1x(2),'go')
plot(r2x(1),r2x(2),'gs')

fprintf(1,'\nLAMBERT TRANSFER:\n\tMin. Delta-v: %f km/s\n\n',...
  minDVTot);

