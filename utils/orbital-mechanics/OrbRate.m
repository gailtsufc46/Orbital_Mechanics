function n = OrbRate( a, varargin )
%
% Compute the mean orbit rate from semi major axis and grav. param. mu
%
% Inputs: 
%   a         Semi major axis 
%   mu        Gravitational parameter (km^3/s^2) 
%
% Outputs:
%   T         Orbital period (sec)
%

% if "mu" not provided, use mu for Earth by default
if( nargin==1 )
  mu = 398600.44;
  n = sqrt(mu./a.^3);
elseif( nargin==3 )
  % if called with 3 inputs, then conform to the SCT usage
  r = varargin{1};
  mu = varargin{2};
  n = sqrt(mu*(2./r - 1./a))./r;
elseif( nargin==2)
  % this is the expected usage, both a and mu provided
  mu = varargin{1};
  n = sqrt(mu./a.^3);
end


