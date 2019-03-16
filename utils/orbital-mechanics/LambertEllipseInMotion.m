function LambertEllipseInMotion( a, c, r1r2 )

% demo
if( nargin < 1 )
  a = 1.25;
  c = 1;
  r1r2 = 3;
end

% pick some initial angle for first drawing
EF = 0;
d = LambertEllipse( a, c, r1r2, EF );

% initialize the plot
handles = LambertEllipsePlot( d );

myslider(handles)

  function myslider(handles)
    
    h = uicontrol(handles.fig,'style','slider','units','pixel','position',[20 20 400 20]);
    try 
        addlistener(h,'ActionEvent',@(hObject, event) makeplot(hObject, event,handles));
    catch
        addlistener(h,'ContinuousValueChange',@(hObject, event) makeplot(hObject, event,handles));
    end
  end
        

  function makeplot(hObject,~,handles)
    
    % get the slider value and translate it into an angle value in radians
    n = get(hObject,'Value');
    EF = n*2*pi;
    
    % compute the new ellipse data
    d = LambertEllipse( a, c, r1r2, EF );
    
    % update the plot data...
    
    set(handles.P1,'xdata',d.x10,'ydata',d.y10);
    set(handles.P2,'xdata',d.x20,'ydata',d.y20);
    set(handles.chord,'xdata',[d.x10 d.x20],'ydata',[d.y10 d.y20]);
    
    % plot the locus of points where the two foci could be
    set(handles.FocusLocus,'xdata',d.xFx,'ydata',d.yFx);
    set(handles.VacantFocusLocus,'xdata',d.xFSx,'ydata',d.yFSx);
    
    
    %% pick an arbitrary main focus
    set(handles.MainFocus,'xdata',d.xF,'ydata',d.yF);
    %plot([0 d.xF],[0 d.yF],'k')
    
    set(handles.VacantFocusLocusP1,'xdata',d.xF1c,'ydata',d.yF1c);
    set(handles.VacantFocusLocusP2,'xdata',d.xF2c,'ydata',d.yF2c);
    
    
    %% plot the exact possible locations of the empty focus
    % its where the circles intersect... note that this intersection is also on
    % the ellipse locus of points that we said the focus had to be on!
    set(handles.VacantFocusPts,'xdata',d.xFS*[1 1],'ydata',[d.yFS1 d.yFS2]);
    
    
    %% plot the ellipses
    set(handles.Ellipse1,'xdata',d.x1,'ydata',d.y1);
    set(handles.Ellipse2,'xdata',d.x2,'ydata',d.y2);
    
    drawnow;
    
  end

end
