function InteractiveTemplate

% A template GUI window with an interactive slider
%

%% Constants
th = linspace(0,2*pi,1e3);

%% Options

%% Default values
a = 3;
b = 2;
c = 0;

%% Init Figure
f = figure('position',[10 10 1000 800],'name','Interactive GUI Template');

ui_handles.a = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .02 .8 .02],...
  'value',a,'min',0,'max',10,'sliderstep',[.01 .05]);

ui_handles.b = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .04 .8 .02],...
  'value',b,'min',0.01,'max',10,'sliderstep',[.01 .05]);

ui_handles.c = uicontrol('parent',f,'style','slider','units','normalized',...
  'position',[.1 .06 .8 .02],...
  'value',c,'min',0.0,'max',2*pi,'sliderstep',[.01 .05]);


% Initialize graphics handles
gfx_handles = InitDisplay;

if( verLessThan('matlab','7.0') )
  actionName = 'ActionEvent';
else
  actionName = 'ContinuousValueChange';
end

fh = @(xx,yy) UpdateDisplay( );
lha = addlistener( ui_handles.a, actionName, fh );
lhb = addlistener( ui_handles.b, actionName, fh );
lhc = addlistener( ui_handles.c, actionName, fh );

  

  function g = InitDisplay

    %% Inertial coordinate axes
    g.xI = line([0 1]*1,[0 0],'color','b','linewidth',2);
    g.yI = line([0 0],[0 1]*1,'color','g','linewidth',2);
    
    if b>a
      sma = b;
      ecc = sqrt(1-a/b);
      phi = pi/2;
    else
      sma = a;
      ecc = sqrt(1-b/a);
      phi = 0;
    end
    
    r = sma*(1-ecc^2) ./ (1+ecc*cos(th));
    xp  = r.*cos(th+phi);
    yp  = r.*sin(th+phi);
    g.path = line(xp,yp,'color','b','linewidth',2);

    
    r  = sma*(1-ecc^2) ./ (1+ecc*cos(c+phi));
    x  = r.*cos(c+phi);
    y  = r.*sin(c+phi);
    g.marker = line(x,y,'color','r','marker','.','markersize',20);
    
    grid on
    zoom on
    axis equal
    axis([-12 12 -12 12])
    
  end

  function UpdateDisplay( )
    
    % Get current values for orbit elements
    a = get(ui_handles.a,'value');
    b = get(ui_handles.b,'value');
    c = get(ui_handles.c,'value');
    
    if b>a      
      sma = b;
      ecc = sqrt(1-a/b);
      phi = pi/2;
    else
      sma = a;
      ecc = sqrt(1-b/a);
      phi = 0;
    end
    
    % Compute new path
    r  = sma*(1-ecc^2) ./ (1+ecc*cos(th));
    xp  = r.*cos(th+phi);
    yp  = r.*sin(th+phi);

    % Compute new point on path
    r  = sma*(1-ecc^2) ./ (1+ecc*cos(c));
    x  = r.*cos(c+phi);
    y  = r.*sin(c+phi);    
            
    % Compute anything else we want to label or show
    
    % Update graphics handles with new data
    set(gfx_handles.path,'xdata', xp,'ydata', yp );
    set(gfx_handles.marker,'xdata',x, 'ydata',y );
    
    axis equal
    axis([-12 12 -12 12])
    
  end


end



