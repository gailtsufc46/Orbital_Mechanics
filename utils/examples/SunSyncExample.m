function SunSyncExample

% Sun sync example
%

%% Constants
mu = 398600.44;
Re = 6378.14;
J2 = 0.00108263;
thx = linspace(0,2*pi,360);

jD0 = J0(2018,4,25);

% orbit parameters 
hp = 800;       % perigee altitude (km)
ha = 800;       % apogee altitude (km)
%inc = 98*pi/180;     % inclination (rad)

W = 0;
w = 0;
th0 = 0;

% compute orbital elements a, e
rp = Re+hp;
ra = Re+ha;
a  = .5*(ra+rp);
e  = abs(ra-rp)/(ra+rp);

WDotDes = 2*pi/365/86400
inc = acos( WDotDes / (-3/2*sqrt(mu)*J2*Re*Re/(1-e^2)^2/a^3.5) );

% compute drift rate
WDot = -3/2*sqrt(mu)*J2*Re*Re/(1-e^2)^2/a^3.5 * cos(inc)
wDot = -3/2*sqrt(mu)*J2*Re*Re/(1-e^2)^2/a^3.5 * (5/2*sin(inc)^2-2);

fprintf(1,'Inclination = %2.1f deg\n',inc*180/pi);

%% Init Figure

% Compute the orbit
thx = linspace(0,2*pi,360);
r = RVFromCOE(a,inc,W,w,e,thx,mu);

% plot the orbit with the Earth
f = figure('position',[10 10 1000 800],'name','Sun-Synchronous Demo');

ui_handles.day = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .02 .8 .02],...
  'value',0,'min',0,'max',365,'sliderstep',[.01 .05]);


% Initialize graphics handles
gfx_handles = InitDisplay;

if( verLessThan('matlab','7.0') )
  actionName = 'ActionEvent';
else
  actionName = 'ContinuousValueChange';
end

fh = @(xx,yy) UpdateDisplay( );
lhday = addlistener( ui_handles.day, actionName, fh );

  

  function g = InitDisplay

    [~,op] = Plot3DOrbit([a,inc,W,w,e,th0],mu,1,[],f);
    g.op = op;
    
    s = SunVector(jD0)*2*Re;    
    g.sun = line([0 s(1)],[0 s(2)],[0 s(3)],'color','y','linewidth',3);
    
    grid on
    zoom on
    axis equal
    
  end

  function UpdateDisplay( )
    
    % Get current values for orbit elements
    day = get(ui_handles.day,'value');
    jD = jD0 + day;
    
    % New values for Right ascension and argument of perigee
    Wt = W + day*86400*WDot;
    wt = w + day*86400*wDot;
    
    % New orbital path
    rt = RVFromCOE(a,inc,Wt,wt,e,thx,mu);
    
    % New sun vector at this day
    s = SunVector(jD)*2*Re;
    
    % Update graphics handles with new data
    set(gfx_handles.op,'xdata', rt(1,:),'ydata', rt(2,:), 'zdata',rt(3,:) );
    set(gfx_handles.sun,'xdata',[0 s(1)], 'ydata',[0 s(2)], 'zdata',[0 s(3)] );
    
    axis equal
    
  end


end








