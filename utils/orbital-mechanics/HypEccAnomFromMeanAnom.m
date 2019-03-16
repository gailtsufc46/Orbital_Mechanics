function F = HypEccAnomFromMeanAnom( Mh, e, tol )
%
% Use Newton Raphson method to solve Keplers equation: e*sinh(F) - F = Mh
%
% Note: The initial guess of F=Mh is only good for relatively small Mh
% values. For larger values of Mh (larger times away from periapsis
% crossing), the value of F should be much smaller than Mh.
%
% Inputs: 
%   Mh    Hyperbolic mean anomaly (rad)
%   e     Eccentricity (0 <= e < 1)
%   tol   Error tolerance. Optional.
%
% Outputs:
%   F     Hyperbolic ecentric anomaly (rad)
%

if nargin<3
  tol = 1e-8;
end

% initial guess 
F = Mh;

dist = inf;
while dist>tol
  
  % compute function value and function derivative value
  f = e*sinh(F)-F-Mh;
  fp = e*cosh(F)-1;

  % distance between current and next iterate
  dist = abs(f/fp);
  
  F = F-f/fp;

end


