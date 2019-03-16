function d = SolarSystem( date, show )

% Load the orbital elements, mass properties, radius and SOI for planets in
% our solar system.
%
% The data in this file was taken from one set of osculating elements
% obtained from the JPL Horizons service with an epoch of 2015-April-06.
%

if( nargin<1 )
  date = datevec(now);
  show = 1;
elseif( nargin<2 )
  show = 1;
end

year = date(1);
month = date(2);
day = date(3);

jd0 = J0(2015,4,6);
dt = (J0(year,month,day) - jd0)*86400; % time in seconds since epoch

%% Constants
sunEarthMass = 332946;
au = 149597871;
G  = 6.673e-11;
mu = 398600.44;
muSun = mu*sunEarthMass;

%% Mercury
%{
Rad = 2440 km
2457118.500000000 = A.D. 2015-Apr-06 00:00:00.0000 (CT)
 EC= 2.056251742522060E-01 QR= 3.075014103054842E-01 IN= 7.004026784982716E+00
 OM= 4.831153048654617E+01 W = 2.917002027982499E+01 Tp=  2457132.325664058793
 N = 4.092339456985709E+00 MA= 3.034206894534501E+02 TA= 2.811216463134211E+02
 A = 3.870986344714715E-01 AD= 4.666958586374588E-01 PR= 8.796924199078170E+01
%}
k = 1;
d(k).name = 'Mercury';
d(k).radius = 2440;
d(k).sma        = 3.870986344714715E-01;
d(k).inc        = 7.004026784982716E+00;
d(k).raan       = 4.831153048654617E+01;
d(k).argPerigee = 2.917002027982499E+01;
d(k).ecc        = 2.056251742522060E-01;
thx             = 2.811216463134211E+02;
d(k).trueAnom   = TrueAnomFromTime( dt, d(k).sma*au, d(k).ecc, muSun, thx*pi/180 )*180/pi;
d(k).earthMass  = .0553;
d(k).mu         = 22032;
d(k).rSOI       = au*d(k).sma*(d(k).earthMass/sunEarthMass)^(2/5);
d(k).el         = [d(k).sma*au, d(k).inc*pi/180, d(k).raan*pi/180, d(k).argPerigee*pi/180, d(k).ecc, d(k).trueAnom*pi/180];

%% Venus
%{
Rad = 6051.8 km
2457118.500000000 = A.D. 2015-Apr-06 00:00:00.0000 (CT)
 EC= 6.759449026890207E-03 QR= 7.184385084418853E-01 IN= 3.394476052282928E+00
 OM= 7.663979781099501E+01 W = 5.467178252963298E+01 Tp=  2457130.893666445743
 N = 1.602145169994792E+00 MA= 3.401435471651922E+02 TA= 3.398783450996168E+02
 A = 7.233278058754326E-01 AD= 7.282171033089799E-01 PR= 2.246987393790104E+02
%}
k = k+1;
d(k).name = 'Venus';
d(k).radius = 6051.8;
d(k).sma        = 7.233278058754326E-01;
d(k).inc        = 3.394476052282928E+00;
d(k).raan       = 7.663979781099501E+01;
d(k).argPerigee = 5.467178252963298E+01;
d(k).ecc        = 6.759449026890207E-03;
thx             = 3.398783450996168E+02;
d(k).trueAnom   = TrueAnomFromTime( dt, d(k).sma*au, d(k).ecc, muSun, thx*pi/180 )*180/pi;
d(k).earthMass  = .815;
d(k).mu         = 324859;
d(k).rSOI       = au*d(k).sma*(d(k).earthMass/sunEarthMass)^(2/5);
d(k).el         = [d(k).sma*au, d(k).inc*pi/180, d(k).raan*pi/180, d(k).argPerigee*pi/180, d(k).ecc, d(k).trueAnom*pi/180];

%% Earth
%{
Rad = 6378 km
2457118.500000000 = A.D. 2015-Apr-06 00:00:00.0000 (CT)
 EC= 1.687449933927703E-02 QR= 9.823213112334106E-01 IN= 6.362225335765567E-04
 OM= 9.209153058977887E+01 W = 8.243575112531003E+00 Tp=  2457023.848318934906
 N = 9.868197179927852E-01 MA= 9.340414521628080E+01 TA= 9.533163022284555E+01
 A = 9.991820073563629E-01 AD= 1.016042703479315E+00 PR= 3.648082759556614E+02
%}
k = k+1;
d(k).name = 'Earth';
d(k).radius = 6378.14;
d(k).sma        = 9.991820073563629E-01;
d(k).inc        = 6.362225335765567E-04;
d(k).raan       = 9.209153058977887E+01;
d(k).argPerigee = 8.243575112531003E+00;
d(k).ecc        = 1.687449933927703E-02;
thx             = 9.533163022284555E+01;
d(k).trueAnom   = TrueAnomFromTime( dt, d(k).sma*au, d(k).ecc, muSun, thx*pi/180 )*180/pi;
d(k).earthMass  = 1;
d(k).mu         = 398600.44;
d(k).rSOI       = au*d(k).sma*(d(k).earthMass/sunEarthMass)^(2/5);
d(k).el         = [d(k).sma*au, d(k).inc*pi/180, d(k).raan*pi/180, d(k).argPerigee*pi/180, d(k).ecc, d(k).trueAnom*pi/180];

%% Mars
%{
Rad = 3389.9 km
2457118.500000000 = A.D. 2015-Apr-06 00:00:00.0000 (CT)
 EC= 9.346452808063631E-02 QR= 1.381206502863846E+00 IN= 1.848408709701127E+00
 OM= 4.951279965997048E+01 W = 2.865383225689425E+02 Tp=  2457003.828511793632
 N = 5.240749793239402E-01 MA= 6.009645781079116E+01 TA= 6.990502762247453E+01
 A = 1.523609991718785E+00 AD= 1.666013480573723E+00 PR= 6.869246085062143E+02
%}
k = k+1;
d(k).name = 'Mars';
d(k).radius = 3389.9;
d(k).sma        = 1.523609991718785E+00;
d(k).inc        = 1.848408709701127E+00;
d(k).raan       = 4.951279965997048E+01;
d(k).argPerigee = 2.865383225689425E+02;
d(k).ecc        = 9.346452808063631E-02;
thx             = 6.990502762247453E+01;
d(k).trueAnom   = TrueAnomFromTime( dt, d(k).sma*au, d(k).ecc, muSun, thx*pi/180 )*180/pi;
d(k).earthMass  = .107;
d(k).mu         = 42828;
d(k).rSOI       = au*d(k).sma*(d(k).earthMass/sunEarthMass)^(2/5);
d(k).el         = [d(k).sma*au, d(k).inc*pi/180, d(k).raan*pi/180, d(k).argPerigee*pi/180, d(k).ecc, d(k).trueAnom*pi/180];

%% Jupiter
%{
Rad = 71492 km
2457118.500000000 = A.D. 2015-Apr-06 00:00:00.0000 (CT)
 EC= 4.898692816046643E-02 QR= 4.946904919016284E+00 IN= 1.303892022825077E+00
 OM= 1.005199597514340E+02 W = 2.737178765816317E+02 Tp=  2455635.815739933401
 N = 8.311715078796857E-02 MA= 1.232364912149122E+02 TA= 1.277743199482675E+02
 A = 5.201721264932293E+00 AD= 5.456537610848303E+00 PR= 4.331236051610554E+03
%}
k = k+1;
d(k).name = 'Jupiter';
d(k).radius = 71492;
d(k).sma        = 5.201721264932293E+00;
d(k).inc        = 1.303892022825077E+00;
d(k).raan       = 1.005199597514340E+02;
d(k).argPerigee = 2.737178765816317E+02;
d(k).ecc        = 4.898692816046643E-02;
thx             = 1.277743199482675E+02;
d(k).trueAnom   = TrueAnomFromTime( dt, d(k).sma*au, d(k).ecc, muSun, thx*pi/180 )*180/pi;
d(k).earthMass  = 318;
d(k).mu         = 126686534;
d(k).rSOI       = au*d(k).sma*(d(k).earthMass/sunEarthMass)^(2/5);
d(k).el         = [d(k).sma*au, d(k).inc*pi/180, d(k).raan*pi/180, d(k).argPerigee*pi/180, d(k).ecc, d(k).trueAnom*pi/180];

%% Saturn
%{
Rad = 60268 km
2457118.500000000 = A.D. 2015-Apr-06 00:00:00.0000 (CT)
 EC= 5.408545737492799E-02 QR= 9.033713158942223E+00 IN= 2.485866924778845E+00
 OM= 1.136359406755999E+02 W = 3.400608608831574E+02 Tp=  2452846.177449937910
 N = 3.339985163983478E-02 MA= 1.426949393296012E+02 TA= 1.462568064350545E+02
 A = 9.550242386455070E+00 AD= 1.006677161396792E+01 PR= 1.077849099097917E+04
%}
k = k+1;
d(k).name = 'Saturn';
d(k).radius = 60268;
d(k).sma        = 9.550242386455070E+00;
d(k).inc        = 2.485866924778845E+00;
d(k).raan       = 1.136359406755999E+02;
d(k).argPerigee = 3.400608608831574E+02;
d(k).ecc        = 5.408545737492799E-02;
thx             = 1.462568064350545E+02;
d(k).trueAnom   = TrueAnomFromTime( dt, d(k).sma*au, d(k).ecc, muSun, thx*pi/180 )*180/pi;
d(k).earthMass  = 95.2;
d(k).mu         = 37931187;
d(k).rSOI       = au*d(k).sma*(d(k).earthMass/sunEarthMass)^(2/5);
d(k).el         = [d(k).sma*au, d(k).inc*pi/180, d(k).raan*pi/180, d(k).argPerigee*pi/180, d(k).ecc, d(k).trueAnom*pi/180];

%% Uranus
%{
Rad = 25559 km
2457118.500000000 = A.D. 2015-Apr-06 00:00:00.0000 (CT)
 EC= 4.919120152409749E-02 QR= 1.821555896848689E+01 IN= 7.712470979255052E-01
 OM= 7.406701967559691E+01 W = 9.656877362142765E+01 Tp=  2470039.705167796463
 N = 1.175410558839414E-02 MA= 2.081227901284173E+02 TA= 2.056035104692815E+02
 A = 1.915796214516051E+01 AD= 2.010036532183414E+01 PR= 3.062759622947914E+04
%}
k = k+1;
d(k).name = 'Uranus';
d(k).radius = 25559;
d(k).sma        = 1.915796214516051E+01;
d(k).inc        = 7.712470979255052E-01;
d(k).raan       = 7.406701967559691E+01;
d(k).argPerigee = 9.656877362142765E+01;
d(k).ecc        = 4.919120152409749E-02;
thx             = 2.056035104692815E+02;
d(k).trueAnom   = TrueAnomFromTime( dt, d(k).sma*au, d(k).ecc, muSun, thx*pi/180 )*180/pi;
d(k).earthMass  = 14.5;
d(k).mu         = 5793939;
d(k).rSOI       = au*d(k).sma*(d(k).earthMass/sunEarthMass)^(2/5);
d(k).el         = [d(k).sma*au, d(k).inc*pi/180, d(k).raan*pi/180, d(k).argPerigee*pi/180, d(k).ecc, d(k).trueAnom*pi/180];

%% Neptune
%{
Rad = 24766 km
2457118.500000000 = A.D. 2015-Apr-06 00:00:00.0000 (CT)
 EC= 8.356298509895565E-03 QR= 2.973091106624503E+01 IN= 1.766339351244111E+00
 OM= 1.317279941885216E+02 W = 2.926004557081222E+02 Tp=  2471454.102766284253
 N = 6.003941020125967E-03 MA= 2.739298865032747E+02 TA= 2.729739367756555E+02
 A = 2.998144497017381E+01 AD= 3.023197887410259E+01 PR= 5.996061566781462E+04
%}
k = k+1;
d(k).name = 'Neptune';
d(k).radius = 24766;
d(k).sma        = 2.998144497017381E+01;
d(k).inc        = 1.766339351244111E+00;
d(k).raan       = 1.317279941885216E+02;
d(k).argPerigee = 2.926004557081222E+02;
d(k).ecc        = 8.356298509895565E-03;
thx             = 2.729739367756555E+02;
d(k).trueAnom   = TrueAnomFromTime( dt, d(k).sma*au, d(k).ecc, muSun, thx*pi/180 )*180/pi;
d(k).earthMass  = 17.2;
d(k).mu         = 6836529;
d(k).rSOI       = au*d(k).sma*(d(k).earthMass/sunEarthMass)^(2/5);
d(k).el         = [d(k).sma*au, d(k).inc*pi/180, d(k).raan*pi/180, d(k).argPerigee*pi/180, d(k).ecc, d(k).trueAnom*pi/180];




%% Plot
if( show )
  PlotSolarSys(d)
end
return;

  function PlotSolarSys(ds)
    
    factor = 1;
    
    f=figure('color','k','name','Solar System');
    
    for i=1:length(ds)
      di = ds(i);
      el = Struct2El(di);
      [a,inc,W,w,e,th0] = OrbitalElements(el);
      thx = linspace(0,2*pi);
      r = RVFromCOE( a,inc,W,w,e,thx, muSun );
      
      % plot the orbit
      figure(f)
      plot3(r(1,:),r(2,:),r(3,:))
      grid on, axis equal, rotate3d on, hold on
      
      r0 = RVFromCOE( a,inc,W,w,e,th0, muSun );
      
      % sphere of influence
      [xs,ys,zs] = sphere(50);
      surf(r0(1)+xs*di.rSOI*factor, r0(2)+ys*di.rSOI*factor, r0(3)+zs*di.rSOI*factor,...
        'edgecolor','none','facealpha',.1,'facecolor','g')
      
      % planet body
      surf(r0(1)+xs*di.radius*factor, r0(2)+ys*di.radius*factor, r0(3)+zs*di.radius*factor,...
        'edgecolor','none','facecolor','r')
      
      plot3(r0(1),r0(2),r0(3),'c.','markersize',20)
      
      u(i) = uicontrol('style','pushbutton','units','normalized',...
        'position',[.05 .9-.1*i .1 .1],'string',di.name,...
        'callback',sprintf('camtarget(%s)',mat2str(r0)));
      
      
    end
    set(gca,'color','k','xcolor','w','ycolor','w','zcolor','w')
    
  end

  function el = Struct2El( dd )
    
    au = 149597871;
    
    el = zeros(1,6);
    el(1) = dd.sma*au;
    el(2) = dd.inc*pi/180;
    el(3) = dd.raan*pi/180;
    el(4) = dd.argPerigee*pi/180;
    el(5) = dd.ecc;
    el(6) = dd.trueAnom*pi/180;
    
  end

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  function j0 = J0(year, month, day)
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %{
  This function computes the Julian day number at 0 UT for any year
  between 1900 and 2100 using Equation 5.48.
 
  j0    - Julian day at 0 hr UT (Universal Time)
  year  - range: 1901 - 2099
  month - range: 1 - 12
  day   - range: 1 - 31
 
  User m-functions required: none
    %}
    % ----------------------------------
    
    j0 = 367*year - fix(7*(year + fix((month + 9)/12))/4) ...
      + fix(275*month/9) + day + 1721013.5;
    
    % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  end %J0

end


