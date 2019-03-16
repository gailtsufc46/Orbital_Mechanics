function InteractiveOrbit

% A GUI to illustrate how an orbit changes as you adjust a few parameters:
%
%   * Semi-major axis
%   * Eccentricity
%   * Inclination
%   * Right Ascension
%   * Argument of Perigee
%   * Mean Anomaly

%% Constants
mu = 398600.44;
Re = 6378.14; % equatorial radius (km)

%% Options
drawEarth = 1;

%% Default orbit values
a = 8000;
e = 0;
inc = 0;
W = 0;
w = 0;
th0 = 0;

%% Init Figure
f = figure('position',[10 10 1000 800],'name','Interactive Orbit GUI');

ui_handles.a = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .02 .8 .02],'value',a,'min',Re,'max',30*Re,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',14,...
'position',[.01 .03 .09 .02],'string','SMA');
ui_handles.aLab = uicontrol('parent',f,'style','text','units','normalized',...
  'position',[.91 .03 .09 .02],'string',num2str(a));

ui_handles.e = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .04 .8 .02],'value',0,'min',0,'max',1,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',14,...
'position',[.01 .05 .09 .02],'string','Ecc.');
ui_handles.eLab = uicontrol('parent',f,'style','text','units','normalized',...
  'position',[.91 .05 .09 .02],'string',num2str(e));

ui_handles.inc = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .06 .8 .02],'value',inc,'min',0,'max',pi,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',14,...
'position',[.01 .07 .09 .02],'string','Inc.');
ui_handles.incLab = uicontrol('parent',f,'style','text','units','normalized',...
  'position',[.91 .07 .09 .02],'string',num2str(inc));

ui_handles.W = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .08 .8 .02],'value',W,'min',0,'max',2*pi,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',14,...
'position',[.01 .09 .09 .02],'string','RAAN');
ui_handles.WLab = uicontrol('parent',f,'style','text','units','normalized',...
  'position',[.91 .09 .09 .02],'string',num2str(W));

ui_handles.w = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .10 .8 .02],'value',w,'min',0,'max',2*pi,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',14,...
'position',[.01 .11 .09 .02],'string','Perigee');
ui_handles.wLab = uicontrol('parent',f,'style','text','units','normalized',...
  'position',[.91 .11 .09 .02],'string',num2str(w));

ui_handles.th0 = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .12 .8 .02],'value',th0,'min',0,'max',2*pi,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',14,...
'position',[.01 .13 .09 .02],'string','True Anom.');
ui_handles.th0Lab = uicontrol('parent',f,'style','text','units','normalized',...
  'position',[.91 .13 .09 .02],'string',num2str(th0));



% Initialize graphics handles
gfx_handles = InitDisplay;

if( verLessThan('matlab','7.0') )
  actionName = 'ActionEvent';
else
  actionName = 'ContinuousValueChange';
end

fh = @(xx,yy) UpdateDisplay( );
lha = addlistener( ui_handles.a, actionName, fh );
lhe = addlistener( ui_handles.e, actionName, fh );
lha = addlistener( ui_handles.inc, actionName, fh );
lhe = addlistener( ui_handles.W, actionName, fh );
lha = addlistener( ui_handles.w, actionName, fh );
lhe = addlistener( ui_handles.th0, actionName, fh );
  

  function g = InitDisplay

    %% Inertial coordinate axes
    g.xI = line([0 1]*2*a,[0 0],[0 0],'color','b','linewidth',2);
    g.yI = line([0 0],[0 1]*2*a,[0 0],'color','g','linewidth',2);
    g.zI = line([0 0],[0 0],[0 1]*2*a,'color','r','linewidth',2);

    %% Draw the Earth (optional)
    if( drawEarth )
      
      image_file = '1024px-Land_ocean_ice_2048.jpg';
      set(gca, 'NextPlot','add', 'Visible','off');
      
      axis equal;
      axis auto;
      axis vis3d;
      hold on;
      
      [x, y, z] = ellipsoid(0, 0, 0, Re, Re, Re, 180);
      g.globe = surf(x, y, -z, 'FaceColor', 'none', 'EdgeColor', 0.5*[1 1 1]);
      
      cdata = imread(image_file);
      set(g.globe, 'FaceColor', 'texturemap', 'CData', cdata, 'FaceAlpha', 1, 'EdgeColor', 'none');
      
      % equatorial plane
      xc = 2*a*cos(0:.01:2*pi); yc = 2*a*sin(0:.01:2*pi);
      g.eq_plane = fill3(xc,yc,xc*0,...
        'c','facealpha',.1);
      
    end
    
    %% Draw the orbit and the orbit plane
    
    % compute the position and velocity from the COE around the orbit
    [th1,th2] = TrueAnomalyRange( e, 1e-1 );
    th_all = linspace(th1,th2,2e3);
    r = RVFromCOE( a,inc,W,w,e,th_all, mu );
    rnow = RVFromCOE( a,inc,W,w,e,th0, mu );
    
    T = OrbPeriod(a,mu);
    t = 0 : T/180 : T*(1-1/180);
    th_t = TrueAnomFromTime(t,a,e,mu,0);
    rt = RVFromCOE(a,inc,W,w,e,th_t,mu);
    
    % plot the orbit
    g.orbit_path = plot3(r(1,:),r(2,:),r(3,:));
    grid on, axis equal, rotate3d on, hold on
    
    g.orbit_now = plot3(rnow(1),rnow(2),rnow(3),'rs','markersize',24,'markerfacecolor','r');    
    g.orbit_pts = plot3(rt(1,:),rt(2,:),rt(3,:),'b.','markersize',16);
    
    % orbital plane
    g.orbit_plane = fill3(r(1,:),r(2,:),r(3,:),...
      'y','facealpha',.5);
    
    
    %% Compute / display apse line and line of nodes for reference
    
    [r0,v0] = RVFromCOE( a, inc, W, w, e, th0, mu );
    [~,~,~,~,~,~,~,lineOfNodes,apseLine] = OrbitalElementsFromRV(r0,v0,mu);
    lineOfNodes = 2*a * lineOfNodes / sqrt(lineOfNodes'*lineOfNodes);
    apseLine    = 2*a * apseLine / sqrt(apseLine'*apseLine);
    angMom = cross(r0,v0)/norm(cross(r0,v0))*a;
    
    g.apse_line = line([0 apseLine(1)],[0 apseLine(2)],[0 apseLine(3)],'color','k','linewidth',2);
    g.apse_line2 = line([0 -apseLine(1)],[0 -apseLine(2)],[0 -apseLine(3)],'color','k','linewidth',2,'linestyle','--');
    g.nodes_line = line([0 lineOfNodes(1)],[0 lineOfNodes(2)],[0 lineOfNodes(3)],'color','m','linewidth',2);
    g.angmom_line = line([0 angMom(1)],[0 angMom(2)],[0 angMom(3)],'color','g','linewidth',2);
    
    g.ra_arc = line([0 0],[0 0],[0 0],'color','b');
    g.inc_arc = line([0 0],[0 0],[0 0],'color','b');
    g.prg_arc = line([0 0],[0 0],[0 0],'color','r');
    
    
    view(120,30)
    camtarget([0 0 0])
    cameratoolbar('show')
    cameratoolbar('setmode','orbit')
    
  end

  function UpdateDisplay( )
    
    % Get current values for orbit elements
    a = get(ui_handles.a,'value');
    e = get(ui_handles.e,'value');
    inc = get(ui_handles.inc,'value');
    W = get(ui_handles.W,'value');
    w = get(ui_handles.w,'value');
    th0 = get(ui_handles.th0,'value');
    
    % update labels
    r2d = 180/pi;
    set(ui_handles.aLab,'string',sprintf('%1.1f',a));
    set(ui_handles.eLab,'string',sprintf('%3.3f',e));
    set(ui_handles.incLab,'string',sprintf('%1.1f',inc*r2d));
    set(ui_handles.WLab,'string',sprintf('%1.1f',W*r2d));
    set(ui_handles.wLab,'string',sprintf('%1.1f',w*r2d));
    set(ui_handles.th0Lab,'string',sprintf('%1.1f',th0*r2d));
    
    % Compute new orbit
    [th1,th2] = TrueAnomalyRange( e, 1e-1 );
    th_all = linspace(th1,th2,5e2);
    r = RVFromCOE( a,inc,W,w,e,th_all, mu );
    
    rnow = RVFromCOE( a,inc,W,w,e,th0, mu );
    
    % Compute faces and vertices of new orbit plane
    plane = patch(r(1,:),r(2,:),r(3,:),'w','visible','off');
    
    T = OrbPeriod(a,mu);
    t = 0 : T/180 : T*(1-1/180);
    th_t = TrueAnomFromTime(t,a,e,mu,0);
    rt = RVFromCOE(a,inc,W,w,e,th_t,mu);
    
    % Compute new apse line, nodes line
    [r0,v0] = RVFromCOE( a, inc, W, w, e, th0, mu );
    [~,~,~,~,~,~,~,lineOfNodes,apseLine] = OrbitalElementsFromRV(r0,v0,mu);
    if( abs(inc)>eps )
      lineOfNodes = 2*a * lineOfNodes / sqrt(lineOfNodes'*lineOfNodes);
    else
      lineOfNodes = [0;0;0];
    end
    if( e>eps )
      apseLine    = 2*a * apseLine / sqrt(apseLine'*apseLine);
    else
      apseLine    = [0;0;0];
    end
    angMom = cross(r0,v0)/norm(cross(r0,v0))*max(Re*2,2*a);
    
    % Compute the inclination, right ascension, and argument of perigee
    % arcs
    % TBD
    ra_arc = GenArc( [1;0;0], lineOfNodes )*Re*2;
    inc_arc = GenArc( [0;0;1], angMom )*Re*2;
    prg_arc = GenArc( lineOfNodes, apseLine )*Re*2;
    
    if(any(imag(inc_arc)))
      1;
    end
    set(gfx_handles.ra_arc,'xdata',ra_arc(1,:),'ydata',ra_arc(2,:),'zdata',ra_arc(3,:));
    set(gfx_handles.inc_arc,'xdata',inc_arc(1,:),'ydata',inc_arc(2,:),'zdata',inc_arc(3,:));
    set(gfx_handles.prg_arc,'xdata',prg_arc(1,:),'ydata',prg_arc(2,:),'zdata',prg_arc(3,:));
    
    % Compute anything else we want to label or show
    
    % Update graphics handles with new data
    set(gfx_handles.orbit_now,'xdata', rnow(1),'ydata', rnow(2),'zdata',rnow(3) );
    set(gfx_handles.orbit_path,'xdata', r(1,:),'ydata', r(2,:),'zdata',r(3,:) );
    set(gfx_handles.orbit_pts,'xdata', rt(1,:),'ydata', rt(2,:),'zdata',rt(3,:) );
    set(gfx_handles.orbit_plane,'faces', plane.Faces,'vertices', plane.Vertices);
    set(gfx_handles.apse_line,'xdata', [0 apseLine(1)] ,'ydata', [0 apseLine(2)],'zdata', [0 apseLine(3)]); 
    set(gfx_handles.apse_line2,'xdata', -[0 apseLine(1)] ,'ydata', -[0 apseLine(2)],'zdata', -[0 apseLine(3)]); 
    set(gfx_handles.nodes_line,'xdata', [0 lineOfNodes(1)],'ydata', [0 lineOfNodes(2)],'zdata', [0 lineOfNodes(3)]); 
    set(gfx_handles.angmom_line,'xdata', [0 angMom(1)],'ydata', [0 angMom(2)],'zdata', [0 angMom(3)]); 
    
  end


end



