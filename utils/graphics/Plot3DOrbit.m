function [f,op] = Plot3DOrbit( el, mu, earth, dT, f )

% Plot an orbit in 3D.
% 
%   Inputs:
%     el      Orbital element vector [a,i,W,w,e,th]
%     mu      Gravitational constant
%     earth   Flag used to show Earth image (1) or not (0).
%     dT      Optional time duration to highlight a segment of the orbit (sec).
%     f       Existing figure handle to add plot to. Optional.
%
%   Outputs:  
%     f       Figure handle for plot
%

if( nargin<1 )
  el = [8000,45,30,0,0,0];
end

if( nargin<2 )
  mu = 398600.44;
end

% orbital elements
[a,inc,W,w,e,th0] = OrbitalElements( el );

%% Draw the Earth (optional)
if( nargin >= 3 && earth )
  
  if( nargin < 5 || ~ishandle(f) )
    f = figure('color', 'k');
  end
  %image_file = 'http://upload.wikimedia.org/wikipedia/commons/thumb/c/cd/Land_ocean_ice_2048.jpg/1024px-Land_ocean_ice_2048.jpg';
  image_file = '1024px-Land_ocean_ice_2048.jpg';
  set(gca, 'NextPlot','add', 'Visible','off');
  
  axis equal;
  axis auto;
  axis vis3d;
  hold on;
  
  Re = 6378.14; % equatorial radius (km)
  [x, y, z] = ellipsoid(0, 0, 0, Re, Re, Re, 180);
  globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
  
  cdata = imread(image_file);
  set(globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');

  % equatorial plane
  xc = 2*a*cos(0:.01:2*pi); yc = 2*a*sin(0:.01:2*pi);
  ep = fill3(xc,yc,xc*0,...
    'c','facealpha',.2);
  
else
  if( nargin < 5 || ~ishandle(f) )
    f = figure('color', 'k');
  end
end

%% Draw the orbit

% compute the position and velocity from the COE around the orbit
[th1,th2] = TrueAnomalyRange( e, 1e-1 );
th_all = linspace(th1,th2);
r = RVFromCOE( a,inc,W,w,e,th_all, mu );

% plot the orbit
figure(f)
plot3(r(1,:),r(2,:),r(3,:))
grid on, axis equal, rotate3d on, hold on

% orbital plane
op = fill3(r(1,:),r(2,:),r(3,:),...
  'y','facealpha',.5);


%% Add vectors for reference
line([0 1]*2*a,[0 0],[0 0],'color','b','linewidth',2)
line([0 0],[0 1]*2*a,[0 0],'color','g','linewidth',2)
line([0 0],[0 0],[0 1]*2*a,'color','r','linewidth',2)

[r0,v0] = RVFromCOE( a, inc, W, w, e, th0, mu );
[~,~,~,~,~,~,~,lineOfNodes,apseLine] = OrbitalElementsFromRV(r0,v0,mu);
lineOfNodes = 2*a * lineOfNodes / sqrt(lineOfNodes'*lineOfNodes);
apseLine    = 2*a * apseLine / sqrt(apseLine'*apseLine);

line([0 apseLine(1)],[0 apseLine(2)],[0 apseLine(3)],'color','y')
line([0 lineOfNodes(1)],[0 lineOfNodes(2)],[0 lineOfNodes(3)],'color','m','linestyle','--')

view(120,30)
camtarget([0 0 0])
cameratoolbar('show')
cameratoolbar('setmode','orbit')

%% Highlight a portion (optional)
if( nargin < 4 || isempty(dT) || dT==0 )
  return
end

% compute the true anomaly values over the given time values "tx"
tx = linspace(0,dT);
thx = TrueAnomFromTime(tx,a,e,mu,th0);

% compute the position and velocity from the COE at each true anomaly "thx"
rx = RVFromCOE( a,inc,W,w,e,thx, mu );

%figure
plot3(rx(1,:),rx(2,:),rx(3,:),'g','linewidth',3)
plot3(rx(1,1),rx(2,1),rx(3,1),'go','linewidth',2,'markersize',14)
plot3(rx(1,end),rx(2,end),rx(3,end),'rs','linewidth',2,'markersize',14)


