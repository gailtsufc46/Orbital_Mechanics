%% Compute the position over time for a hyperbolic trajectory

% load the data from problem_3_18
problem_3_18

% define a time vector from the starting point to a mirror point on the
% other side of Earth, sample every 10 minutes
tx = t2 : 60 : -t2;

% compute the hyperbolic mean anomaly at each time point.
Mhx = (mu^2/h^3)*(e^2-1)^1.5 * tx;

% compute the hyperbolic eccentric anomaly at each point
nx = length(Mhx);
Fx = zeros(1,nx);
tol = 1e-8;
for i=1:nx
  Fx(i) = HypEccAnomFromMeanAnom( Mhx(i), e, tol );
end

% compute the true anomaly at each point
thx = TrueAnomFromHypEccAnom( Fx, e );

% compute the radius at each point
rx = h^2/mu ./ (1+e*cos(thx));
% -or, equivalently-
rx2 = a*(e*cosh(Fx)-1);

% coordinates in the perifocal frame
x = rx.*cos(thx);
y = rx.*sin(thx);

% plot the trajectory
figure, plot(x,y,'-s'), axis equal, grid on, hold on
set(gca,'fontsize',18)

% add Earth
ang = linspace(0,2*pi);
fill(Re*cos(ang),Re*sin(ang),'b')

%% Animate the motion along the trajectory

% compute the speed at each point
energy = v0^2/2 - mu/r0;
vx = sqrt( 2*(energy+mu./rx) );

sat = line('xdata',x(1),'ydata',y(1),'marker','o','color','r',...
  'linewidth',3,'markersize',20);
%set( sat, 'erasemode','xor' );

speedText = text(x(1),y(1),sprintf('V = %2.1f km/s',vx(1)));
%set(speedText,'erasemode','xor','fontsize',14);

speedFactor = 5e3;
waitTime = (tx(2)-tx(1))/speedFactor;

figure(gcf)
fig = msgbox('Click to Start');
uiwait(fig)
for i=1:nx
  set(sat,'xdata',x(i),'ydata',y(i))
  set(speedText,'position',[x(i),y(i),0],...
    'string',sprintf('  V = %2.1f km/s',vx(i)),'fontsize',24)
  pause(waitTime)
end




