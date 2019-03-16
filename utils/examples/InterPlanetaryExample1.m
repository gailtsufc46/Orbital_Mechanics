%% Example of an interplanetary transfer
%
%   Earth to Mars transfer. 
%   
%   Planetary Departure
%     Start in LEO parking orbit
%     Delta-V to get onto hyperbolic orbit.
%     Hyperbolic orbit out to SOI
% 
%   Hohmann transfer coast from Earth SOI to Mars SOI
% 
%   Planetary Rendezvous
%     Enter Mars SOI on hyperbolic orbit
%     Delta-V at perigee to circularize
%     End in Mars parking orbit

%% Constants
au    = 149597871;
rE    = 149.6e6; %au;
rM    = 227.9e6; %1.5*au;
muE   = 398600.44;
muM   = 42828;
muS   = 1.327e11; %132712440018;
radE  = 6378.14;
radM  = 3389.90;
rSOIE = 923890.49;
rSOIM = 576239.75;

%% Hohmann Transfer
dVHE  = sqrt(muS/rE)*( sqrt(2*rM/(rE+rM)) - 1 );
dVHM  = sqrt(muS/rM)*( 1 - sqrt(2*rE/(rE+rM)) );


%% Planetary Departure at Earth
rP1   = radE+300;                       % parking orbit radius and perigee of hyperbolic transfer
vC1   = sqrt(muE/rP1);                  % circular orbit velocity
vInf1 = dVHE;                           % hyperbolic excess velocity is same as delta-v for hohmann
h1    = rP1*sqrt(vInf1^2 + 2*muE/rP1);  % angular momentum of hyperbola
e1    = 1+rP1*vInf1^2 / muE;            % eccentricity of hyperbola
vP1   = h1/rP1;                         % velocity of hyperobola at perigee 
dV1   = vP1-vC1;                        % required delta-v from parking orbit
beta1 = acos(1/e1);                     % angle between apse line and asymptote line for hyperbola
alpha1  = pi/2-beta1;                   % angel between apse line and line to sun
sma1    = h1^2/muE/(1-e1^2);            % semi major axis of hyperbola
delta1  = (rP1+abs(sma1))*sin(beta1);   % perp. distance between asymptote and line parallel to asymptote that goes through Earth center


%% Planetary Arrival at Mars
vInf2 = dVHM;
T2    = 7*3600;  % period of final orbit at Mars
aF    = (T2*sqrt(muM)/2/pi)^(2/3);
eF    = 2*muM/aF/vInf2^2-1; % found by solving for e using optimal periapse formula and rp as f(a,e)
hF    = sqrt(aF*muM*(1-eF^2));
rP2   = aF*(1-eF);
h2    = rP2*sqrt(vInf2^2 + 2*muM/rP2);
vP2   = h2/rP2;
dV2   = vP2-hF/rP2;
delta2 = rP2*sqrt(2/(1-eF));
beta2  = acos( 1/(1+rP2*vInf2^2/muM) );



