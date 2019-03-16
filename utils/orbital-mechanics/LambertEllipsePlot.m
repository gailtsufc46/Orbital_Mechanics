function handles = LambertEllipsePlot( d )

% Plot the ellipses that pass through two points P1 and P2, based on
% Lambert's theorem

% demo
if( nargin < 1 )
  a = 1.25;
  c = 1;
  r1r2 = 2;
  EF = 1.63;
  d = LambertEllipse(a,c,r1r2,EF);
end

% plot points P1 and P2, and the chord between them
handles.fig = figure('name','Lambert Ellipse Plot','position',[20 20 800 1200]);
handles.P1 = plot(d.x10,d.y10,'b.','markersize',40,'linewidth',3); 
hold on, axis equal
handles.P2 = plot(d.x20,d.y20,'r.','markersize',40,'linewidth',3);
handles.chord = plot([d.x10 d.x20],[d.y10 d.y20],'k--');

% plot the locus of points where the two foci could be
handles.FocusLocus = plot(d.xFx,d.yFx,'k');
handles.VacantFocusLocus = plot(d.xFSx,d.yFSx,'m');


%% pick an arbitrary main focus 
handles.MainFocus = plot(d.xF,d.yF,'ks','markersize',20,'linewidth',3);
%plot([0 d.xF],[0 d.yF],'k')

handles.VacantFocusLocusP1 = plot(d.xF1c,d.yF1c,'b');
handles.VacantFocusLocusP2 = plot(d.xF2c,d.yF2c,'r');


%% plot the exact possible locations of the empty focus
% its where the circles intersect... note that this intersection is also on
% the ellipse locus of points that we said the focus had to be on!
handles.VacantFocusPts = plot(d.xFS*[1 1],[d.yFS1 d.yFS2],'k*');


%% plot the ellipses
handles.Ellipse1 = plot(d.x1,d.y1,'g--','linewidth',2);
handles.Ellipse2 = plot(d.x2,d.y2,'c:','linewidth',2);



