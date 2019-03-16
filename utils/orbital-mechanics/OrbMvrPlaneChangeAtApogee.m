function d = OrbMvrPlaneChangeAtApogee( v1, v2, delta, method )

%
% Compute an orbit maneuver for a plane change at apogee
%
% Inputs: 
%   v1      Initial speed
%   v2      Final speed
%   delta   Plane rotation angle (rad)
%   method  'simultaneous' OR 'plane-speed' OR 'speed-plane'
%
% Outputs:
%   d     Data structure with fields:
%         .dVspeed
%         .dVplane
%         .dVTot
%

if( nargin < 4 )
  method = 'simultaneous';
  warning('No method entered. Using "simultaneous".');
end

switch lower(method)
  case {1,'simultaneous','i'}
    
    d.dVTot = sqrt( (v2-v1)^2 + 4*v1*v2*sin(delta/2)^2 );
    
  case {2,'plane-speed','ii'}
    
    d.dVspeed = abs(v2-v1);
    d.dVplane = 2*v1*sin(delta/2);
    d.dVTot = d.dVspeed + d.dVplane;
    
  case {3,'speed-plane','iii'}
    
    d.dVspeed = abs(v2-v1);
    d.dVplane = 2*v2*sin(delta/2);
    d.dVTot = d.dVspeed + d.dVplane;
    
  otherwise
    error('Method not recognized. Must use ''simultaneous'', ''speed-plane'', or ''plane-speed''.');
end
