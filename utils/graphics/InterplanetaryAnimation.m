function d = InterplanetaryAnimation( planets, sc, jD0 )

% Animate an interplanetary mission
%
%

%% Input checking and default values
if( nargin<1 )
  d = InterplanetaryAnimation( [2 3 4] );
  return
end

if( nargin<2 )
  sc = [];
end

if( nargin<3 )
  jD0 = now + 1721058.5;
end

%% Constants
mu = 1.327e11;

%% Time
date = datevec(jD0-1721058.5);
yr = date(1);
mo = date(2);
day = date(3);

if( isempty(sc) )
  time = (0 : 1 : 365*12)*86400;
else
  
  % get all of the delta-v dates
  dVJD = [sc.dVJD];
  
  % Define the time vector to go all the way to the last delta-v date, plus
  % one Earth year, with intervals of 1 day.
  time = (0 : 1 : (max(dVJD)-jD0+365))*86400;
  
  % Now insert all of the delta-v times, sort, and remove duplicates
  dVTimes = (dVJD-jD0)*86400;
  time = unique(sort([time, dVTimes]));
  
end

%% Initialize the Figure
f1=figure;
plot(0,0,'k.','markersize',30), hold on, axis equal, grid on
tt = title('','fontsize',18);

%% Plot the orbits of the planets and their initial positions

for i=1:length(planets)
  [r0,v0,~,coe] = PlanetData( planets(i), yr, mo, day, 0, 0, 0, mu );
  T = 2*pi/sqrt(mu/coe(1)^3);
  thx = linspace(0,2*pi);
  rx = RVFromCOE(coe(1),coe(2),coe(3),coe(4),coe(5),coe(6)+thx,mu);
  plot3(rx(1,:),rx(2,:),rx(3,:)), hold on
  hp(i) = plot3(rx(1,1),rx(2,1),rx(3,1),'.','markersize',24);
end
axis equal

%% Plot the initial positions of the planets and the spacecraft
for i=1:length(sc)
  hs(i) = plot3(sc(i).r0(1),sc(i).r0(2),sc(i).r0(3),'.','markersize',24);
end

ylim=get(gca,'ylim'); xlim=get(gca,'xlim');
dvText = text(xlim(1)+0.9*diff(xlim),ylim(1)+0.9*diff(ylim),'');

%% Add a slider
timeSlider = uicontrol('parent',f1,'style','slider','units','normalized',...
  'position',[.1 .05 .8 .02],'min',0,'max',time(end),'value',0,'sliderstep',[.0001 .01]);

if( verLessThan('matlab','7.0') )
  actionName = 'ActionEvent';
else
  actionName = 'ContinuousValueChange';
end

fh = @(xx,yy) UpdateDisplay( );
lha = addlistener( timeSlider, actionName, fh );

%% Compute and store the trajectory data

% first we need a placeholder for planet position and spacecraft position data
rp = cell(1,length(planets));
rs = cell(1,length(sc));
d = ComputeTrajectories;


%% Compute Trajectories
  function d = ComputeTrajectories
    
    for i=1:length(planets)
      [r0p,v0p] = PlanetData(planets(i),yr,mo,day,0,0,0,mu);
      rp{i} = RVAtTFromR0V0(r0p,v0p,time,mu);
    end
    
    lastTime = -1;
    for i=1:length(sc)
      r0s = sc(i).r0;
      v0s = sc(i).v0;
      for j=1:length(sc(i).dVJD)
        nextTime = (sc(i).dVJD(j)-jD0)*86400;
        tx = time(time>lastTime & time<=nextTime);
        [rsx,vsx] = RVAtTFromR0V0(r0s,v0s,tx-tx(1),mu);
        r0s = rsx(:,end);
        v0s = vsx(:,end) + sc(i).dV(:,j);
        rs{i} = [rs{i},rsx];
        lastTime = nextTime;
      end
      nextTime = time(end);
      tx = time(time>lastTime & time<=nextTime);
      [rsx,vsx] = RVAtTFromR0V0(r0s,v0s,tx-tx(1),mu);
      rs{i} = [rs{i},rsx];
      
    end
    
    d.planets = planets;
    d.sc = sc;
    d.rp = rp;
    d.rs = rs;
    d.t  = time;
    d.jD0 = jD0;
    s.mu = mu;
    
  end


%% Update the display
  function UpdateDisplay
    
    tk = get(timeSlider,'value');
    set(tt,'string',datestr(jD0+tk/86400-1721058.5,'dd-mmm-yyyy'));
    
    % Compute the new position of each planet
    for i=1:length(planets)
      rpt = interp1(time,rp{i}',tk);
      set(hp(i),'xdata',rpt(1),'ydata',rpt(2),'zdata',rpt(3));
    end
    
    % Compute the new position of each spacecraft
    for i=1:length(sc)
      rst = interp1(time,rs{i}',tk);
      set(hs(i),'xdata',rst(1),'ydata',rst(2),'zdata',rst(3));
    end
    
    % denote times when delta-v's occur
    if( ~isempty(sc) )
      dvjd = [sc.dVJD];
      k = find( abs(dvjd-(jD0+tk/86400))<=5 );
    else
      k=[];
    end
    if( ~isempty(k) )
      set(gca,'color',[1 .5 .5])
      set(dvText,'string',sprintf('Delta-V %d',k));
    else
      set(gca,'color','w')
      set(dvText,'string','');
    end
  end

end







