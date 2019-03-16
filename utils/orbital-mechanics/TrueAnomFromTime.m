function th = TrueAnomFromTime( t, a, e, mu, th0 )
%
% Compute the true anomaly given the time elapsed since th0
% Note: If th0 is not supplied, we set th0 = 0, so it is equivalent to
%       the time since periapse crossing.
%
% USAGE:
%   th = TrueAnomFromTime( t, a, e, mu, th0 )
%
% Inputs: 
%   t         Time measured since periapsis crossing (seconds) 
%   a         Semi major axis 
%   e         Eccentricity 
%   mu        Gravitational parameter (km^3/s^2) 
%
% Outputs:
%   th        True anomaly (rad)
%

if( nargin<5 )
  th0 = 0;
end

n = OrbRate(a,mu);
th = zeros(size(t));

% compute the true anomaly values "th" over all time points "t"
if( e<1 )
  
  % first compute the initial mean anomaly 
  E0 = EccAnomFromTrueAnom(th0,e);
  M0 = MeanAnomFromEccAnomE(E0,e);

  % now compute the mean anomaly over time
  M = M0+n*t;
  
  % finally convert back to true anomaly
  for i=1:length(t)
    E = EccAnomFromMeanAnom(M(i),e);
    th(i) = TrueAnomFromEccAnom(E,e);
  end
  
elseif( e==1 )
  % TBD
else
  
  % first compute the initial hyperbolic mean anomaly
  F0 = HypEccAnomFromTrueAnom(th0,e);
  Mh0 = HypMeanAnomFromHypEccAnom(F0,e);
  
  % next compute the hyperbolic mean anomaly over time
  h = sqrt(a*mu*(1-e^2));
  Mh = Mh0+(mu^2/h^3)*(e^2-1)^1.5*t;
  
  % finally convert back to true anomaly
  for i=1:length(t)
    F = HypEccAnomFromMeanAnom(Mh(i),e);
    th(i) = TrueAnomFromHypEccAnom(F,e);
  end
  
end


