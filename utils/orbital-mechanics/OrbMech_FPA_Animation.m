function OrbMech_FPA_Animation

% Show how the flight path angle varies across the orbit
%
% For use with AEM 4301 "Orbital Mechanics"
% Author: Joseph Mueller


%% Constants
th = linspace(0,2*pi,1e3);

%% Options

%% Default values
a = 3;
e = .25;
ta = pi/4;

%% Init Figure
f = figure('position',[10 10 1000 800],'name','Interactive GUI Template');

ui_handles.a = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .01 .8 .03],...
  'value',a,'min',1,'max',10,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',16,...
  'position',[.01 .02 .09 .03],'string','SMA');
ui_handles.a_lab = uicontrol('parent',f,'style','text','units','normalized',...
  'fontsize',16,'position',[.91 .02 .09 .03],'string',sprintf('%2.1f',a));

ui_handles.e = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .04 .8 .03],...
  'value',e,'min',0.05,'max',.95,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',16,...
  'position',[.01 .05 .09 .03],'string','Ecc');
ui_handles.e_lab = uicontrol('parent',f,'style','text','units','normalized',...
  'fontsize',16,'position',[.91 .05 .09 .03],'string',sprintf('%4.4f',e));

ui_handles.ta = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .07 .8 .03],...
  'value',ta,'min',0.0,'max',4*pi,'sliderstep',[.01 .05]);
uicontrol('parent',f,'style','text','units','normalized','fontsize',16,...
  'position',[.01 .08 .09 .03],'string','T.A.');
ui_handles.ta_lab = uicontrol('parent',f,'style','text','units','normalized',...
  'fontsize',16,'position',[.91 .08 .09 .03],'string',sprintf('%2.1f deg',ta*180/pi));


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
lhta = addlistener( ui_handles.ta, actionName, fh );

  

  function g = InitDisplay

    g.axes = axes('parent',f,'units','normalized','position',[0.13 0.15 0.775 0.8]);
    
    %% Inertial coordinate axes
    g.xI = line([0 1]*1,[0 0],'color','b','linewidth',2);
    g.yI = line([0 0],[0 1]*1,'color','g','linewidth',2);
    
    hold on
    
    rp = a*(1-e^2) ./ (1+e*cos(th));
    xp  = rp.*cos(th);
    yp  = rp.*sin(th);
    g.path = line(xp,yp,'color','k','linewidth',3);

    
    r  = a*(1-e^2) ./ (1+e*cos(ta));
    x  = r.*cos(ta);
    y  = r.*sin(ta);
    g.marker = line(x,y,'color','r','marker','.','markersize',20);
    
    % Compute vr and vt
    h = sqrt(a*(1-e^2));
    vr = e/h*sin(ta);
    vt = h/r;
    vx = vr*cos(ta)-vt*sin(ta);
    vy = vr*sin(ta)+vt*cos(ta);
    vv = [vx;vy];
    g.vel = line([x, x+5*vx],[y, y+5*vy],'color','m','linewidth',1.5);
    
    dirLen = a*.7;
    
    g.rad   = line('xdata',[x, x+dirLen*cos(ta)], 'ydata',[y, y+dirLen*sin(ta)],'color','k','linestyle','--');    
    g.trans = line('xdata',[x, x-dirLen*sin(ta)], 'ydata',[y, y+dirLen*cos(ta)],'color','k','linestyle','--');    
    g.r = line('xdata',[0,x], 'ydata',[0,y],'color','k','linestyle','--');    
    
    % Gamma label in title
    gam = atan2d(vr,vt);
    g.gammaText = title(sprintf('\\gamma = %2.2f deg',gam),'fontsize',18);
 
    % Angular sweep with patch
    n = 100;
    fac = linspace(0,1,n);
    xs = zeros(1,n);
    ys = xs;
    vtv = vt*[-sin(ta);cos(ta)];
    for i=1:n
      vs=vtv*(1-fac(i))+vv*fac(i);
      vs=vs/norm(vs)*a/2;
      xs(i)=x+vs(1);
      ys(i)=y+vs(2);
    end
    if( gam>0 )
      col = 'g';
    else
      col = 'r';
    end
    g.gammaSweep = fill([x,xs,x],[y,ys,y],col,'facealpha',.5,'edgecolor','none');

    
    grid on
    zoom on
    axis equal
    axis([min(xp)-3, max(xp)+3, min(yp)-3, max(yp)+3])
    
  end

  function UpdateDisplay( )
    
    % Get current values for orbit elements
    a = get(ui_handles.a,'value');
    e = get(ui_handles.e,'value');
    ta = get(ui_handles.ta,'value');
        
    % Compute new path
    rp  = a*(1-e^2) ./ (1+e*cos(th));
    xp  = rp.*cos(th);
    yp  = rp.*sin(th);

    % Compute new point on path
    r  = a*(1-e^2) ./ (1+e*cos(ta));
    x  = r.*cos(ta);
    y  = r.*sin(ta);    
            
    % Compute vr and vt and vx and vy
    h = sqrt(a*(1-e^2));
    vr = e/h*sin(ta);
    vt = h/r;
    vx = vr*cos(ta)-vt*sin(ta);
    vy = vr*sin(ta)+vt*cos(ta);
    vv = [vx;vy];
        
    gam = atan2d(vr,vt);
    
    dirLen = a*.7;    
    
    % Angular sweep with patch
    n = 100;
    fac = linspace(0,1,n);
    xs = zeros(1,n);
    ys = xs;
    vtv = vt*[-sin(ta);cos(ta)];
    for i=1:n
      vs=vtv*(1-fac(i))+vv*fac(i);
      vs=vs/norm(vs)*a/2;
      xs(i)=x+vs(1);
      ys(i)=y+vs(2);
    end
    dummy = patch([x xs x],[y ys y],'w','visible','off');
    %figure, patch([x xs x],[y ys y],'r');
    
    % Update graphics handles with new data
    set(gfx_handles.path,'xdata', xp,'ydata', yp );
    set(gfx_handles.marker,'xdata',x, 'ydata',y );
    set(gfx_handles.vel,'xdata',[x, x+5*vx], 'ydata',[y, y+5*vy]);
    set(gfx_handles.gammaText,'string',sprintf('\\gamma = %2.2f deg',gam));
    
    set(gfx_handles.rad,'xdata',[x, x+dirLen*cos(ta)], 'ydata',[y, y+dirLen*sin(ta)]);
    set(gfx_handles.trans,'xdata',[x, x-dirLen*sin(ta)], 'ydata',[y, y+dirLen*cos(ta)]);
    set(gfx_handles.r,'xdata',[0, x], 'ydata',[0, y]);
    set(gfx_handles.gammaSweep,'faces',dummy.Faces,'vertices',dummy.Vertices);
    
    if( gam>0 )
      set(gfx_handles.gammaSweep,'facecolor','g')
    else
      set(gfx_handles.gammaSweep,'facecolor','r');
    end
    
    set(ui_handles.a_lab,'string',sprintf('%2.1f',a));
    set(ui_handles.e_lab,'string',sprintf('%4.4f',e));
    set(ui_handles.ta_lab,'string',sprintf('%2.1f deg',ta*180/pi));
    
    axis equal
    axis([min(xp)-3, max(xp)+3, min(yp)-3, max(yp)+3])
    
  end


end



