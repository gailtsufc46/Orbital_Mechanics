function [tW,tT,tSyn] = HohmannTransferAnimation( r1, r2, th10, th20, mu )

% Animate a Hohmann Transfer Rendezvous
%
%


if( nargin<1 )
  mu = 1.32712438e11;
  s = SolarSystem(0);
  k1 = 4;
  k2 = 3;
  r1 = s(k1).el(1);
  r2 = s(k2).el(1);
  th10 = (s(k1).raan + s(k1).argPerigee + s(k1).trueAnom)*pi/180;
  th20 = (s(k2).raan + s(k2).argPerigee + s(k2).trueAnom)*pi/180;
  HohmannTransferAnimation( r1, r2, th10, th20, mu );
  return
end

%% Compute the wait time, maneuver time, and all locations
v1 = sqrt( mu/r1 );
v2 = sqrt( mu/r2 );

n1 = sqrt( mu/r1^3 );
n2 = sqrt( mu/r2^3 );

hT = sqrt( 2*mu * r1*r2 / (r1+r2) );
tT = pi/sqrt(mu)*((r1+r2)/2)^(3/2);

tSyn = 2*pi/abs(n2-n1);
tW = (th20 - th10 + tT*n2 + pi)/(n1-n2);
while( tW < 0 )
  tW = tW + tSyn;
end

% TA1 at start of transfer
th1s = th10 + tW*n1;

% TA2 at start of transfer
th2s = th20 + tW*n2;

% TA1 at end of transfer
th1f = th1s + tT*n1;

% TA2 at end of transfer
th2f = th2s + tT*n2;

%% Plot the orbits and all locations
f1=figure;
plot(0,0,'k.','markersize',30), hold on, axis equal, grid on
tt = title('');

thx = linspace(0,2*pi,300);
plot(r1*cos(thx),r1*sin(thx),'b');
plot(r2*cos(thx),r2*sin(thx),'r');

dT = tT/1e4;

twx = 0:dT:tW;
tmx = tW : dT : (tW+tT);
tpx = (tW+tT+.01) : dT : (tW+tT+2*pi/n2+.01);
t = [twx,tmx,tpx];

aT = .5*(r1+r2);
eT = abs(r2-r1)/(r2+r1);
if( r1 < aT ) % inner to outer (transfer starts at perigee)
  thT = TrueAnomFromTime(tmx-tW,aT,eT,mu,0);
  M = [cos(th1s), -sin(th1s), 0; sin(th1s), cos(th1s), 0; 0 0 1];
else
  % outer to inner (transfer starts at apogee)
  thT = TrueAnomFromTime(tmx-tW,aT,eT,mu,pi);
  M = [cos(th1s+pi), -sin(th1s+pi), 0; sin(th1s+pi), cos(th1s+pi), 0; 0 0 1];
end
[rT,xT,yT] = PerifocalOrbit( aT, eT, thT );
rrT = M*[xT; yT; 0*rT];

plot(rrT(1,:),rrT(2,:),'k')
h1=plot(r1*cos(th10),r1*sin(th10),'b.','markersize',24);
h2=plot(r2*cos(th20),r2*sin(th20),'r.','markersize',24);
hs=plot(r1*cos(th10),r1*sin(th10),'mo','linewidth',2,'markersize',16);
plot(r1*cos(th1s),r1*sin(th1s),'b*')
plot(r2*cos(th2s),r2*sin(th2s),'r*')
plot(r1*cos(th1f),r1*sin(th1f),'b^')
plot(r2*cos(th2f),r2*sin(th2f),'r^')
plot(r1*cos(th1s),r1*sin(th1s),'ko','markersize',18)
plot(r2*cos(th1s+pi),r2*sin(th1s+pi),'ks','markersize',18)

%% Plot a time history of x and y for spacecraft 1 and planet 2
x1 = r1*cos(th10+n1*t);
y1 = r1*sin(th10+n1*t);
x2 = r2*cos(th20+n2*t);
y2 = r2*sin(th20+n2*t);
xs = [r1*cos(th10+n1*twx), rrT(1,:), r2*cos(th1s+pi+n2*(tpx-tpx(1)))];
ys = [r1*sin(th10+n1*twx), rrT(2,:), r2*sin(th1s+pi+n2*(tpx-tpx(1)))];

figure
subplot(211)
plot(t,r2*cos(th20+n2*t),'k'), hold on, grid on
plot(twx,r1*cos(th10+n1*twx),tmx,rrT(1,:),tpx,r2*cos(th1s+pi+n2*(tpx-tpx(1))),'--')
ylabel('x (km)')
subplot(212)
plot(t,r2*sin(th20+n2*t),'k'), hold on, grid on
plot(twx,r1*sin(th10+n1*twx),tmx,rrT(2,:),tpx,r2*sin(th1s+pi+n2*(tpx-tpx(1))),'--')
ylabel('y (km)')
xlabel('Time (sec)')


%% Add a slider
timeSlider = uicontrol('parent',f1,'style','slider','units','normalized',...
  'position',[.1 .05 .8 .02],'min',0,'max',t(end),'value',0,'sliderstep',[.0001 .01]);

if( verLessThan('matlab','7.0') )
  actionName = 'ActionEvent';
else
  actionName = 'ContinuousValueChange';
end

fh = @(xx,yy) UpdateDisplay( );
lha = addlistener( timeSlider, actionName, fh );

%% Update the display
  function UpdateDisplay
    
    tk = get(timeSlider,'value');
    
    set(tt,'string',sprintf('Day %d',round(tk/86400)))
    
    x1t=interp1(t,x1,tk);
    y1t=interp1(t,y1,tk);
    x2t=interp1(t,x2,tk);
    y2t=interp1(t,y2,tk);
    xst=interp1(t,xs,tk);
    yst=interp1(t,ys,tk);
    
    set(h1,'xdata',x1t,'ydata',y1t);
    set(h2,'xdata',x2t,'ydata',y2t);
    set(hs,'xdata',xst,'ydata',yst);
    
  end

end






