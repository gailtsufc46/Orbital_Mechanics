% Problem 4.2
% From "Orbital Mechanics for Engineering Students" 3rd Ed., Curtis.

%% given
Re = 6378;    % km
mu = 398600;  % km^3/s^2
h0 = 500;     % km
v0 = 10;      % km/s, in the K direction
RA = 300*pi/180;  % rad
dec = -20*pi/180; % rad
dT = 1800;        % sec

%% find RA and dec 30 min later - using the method I described in class

% point 0 corresponds to the initial time, with given pos. and vel. above
% point 1 corresponds to the position 30 min later ...

% Compute position magnitude, position vector, velocity vector in ECI
r0 = Re+h0;
rvec = r0*[cos(dec)*cos(RA); cos(dec)*sin(RA); sin(dec)];
vvec = [0; 0; v0];

% Now we can compute several properties of the orbit ...

vr = vvec'*rvec/r0;               % radial velocity
vt = sqrt(v0^2-vr^2);             % transverse velocity
h  = vt*r0;                       % angular momentum (specific)
theta0 = atan2( vr*h, vt*h-mu );  % initial true anomaly
e = vr*h/mu/sin(theta0);          % eccentricity
a = h^2 / mu / (1-e^2);           % semi major axis
n = sqrt(mu/a^3);                 % mean orbit rate

% what is the initial mean anomaly?
E0 = EccAnomFromTrueAnom(theta0,e);
M0 = MeanAnomFromEccAnomE( E0, e );

% what is the new mean anomaly after dT? This is at t1 = t0 + dT
dM = n*dT;
M1 = M0 + dM;

% now compute the true anomaly at time t1
E1 = EccAnomFromMeanAnom(M1,e);
theta1 = TrueAnomFromEccAnom(E1,e);

% compute absolute value of change in true anomaly
deltaTA = theta1-theta0;
if( deltaTA<0 )
  deltaTA = 2*pi+deltaTA;
end

% since the orbital plane is orthogonal to the equatorial plane...
%{
 * declination changes uniformly with change in true anomaly
      - ARG = dec0 + abs(theta1-theta0)
      - dec1 = ARG      if cos(ARG)>0
      - dec1 = 180-ARG  if cos(ARG)<0
 * right ascension flips back and forth between 300 and 120 
      - RA1 = 300       if cos(ARG)>0
      - RA1 = 120       if cos(ARG)<0
%}

% new right ascension and declination should be:
ARG = dec+deltaTA;
if( cos(ARG)<0 )
  dec1 = pi-ARG;
  RA1 = 120*pi/180;
else
  dec1 = ARG;
  RA1 = 300*pi/180;
end

%% find RA and dec 30 min later - use Alg. 3.4 and Alg. 4.1
[rvec1b,vvec1b]=rv_from_r0v0(rvec,vvec,dT); % from App. D.
[RA1b,dec1b] = ra_and_dec_from_r(rvec1b);

fprintf(1,'According to our method...: RA = %2.2f deg, dec = %2.2f deg\n',RA1*180/pi,dec1*180/pi)
fprintf(1,'According to Algorithm 3.4: RA = %2.2f deg, dec = %2.2f deg\n',RA1b,dec1b)

%% plot the orbit and the initial / final points on it

% position at new point dT later
r1mag = a*(1-e^2)/(1+e*cos(theta1));
rvec1 = r1mag*[cos(dec1)*cos(RA1); cos(dec1)*sin(RA1); sin(dec1)];

nn=1e4; 
th = linspace(0,2*pi,nn); 
rr = zeros(3,nn); 
vv=zeros(3,nn); 
for i=1:nn, 
  [rr(:,i),vv(:,i)] = rv_from_r0v0_ta(rvec, vvec, th(i)*180, mu); 
end
figure, 
plot3(rr(1,:),rr(2,:),rr(3,:)), grid on, axis equal
hold on
plot3(0,0,0,'b.','markersize',24)
plot3(rvec(1),rvec(2),rvec(3),'go','markersize',20,'linewidth',2)
plot3(rvec1(1),rvec1(2),rvec1(3),'rs','markersize',20,'linewidth',2)
plot3(rvec1b(1),rvec1b(2),rvec1b(3),'m*')

% plot Earth
[x,y,z] = sphere(250);
hold on
surf(x*Re,y*Re,z*Re,'facecolor',[0 0 .8],'edgecolor','none','facealpha',.1)



