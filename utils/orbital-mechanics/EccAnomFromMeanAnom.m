function E = EccAnomFromMeanAnom( Me, e, tol )
%
% Use the Newton Raphson method to solve Keplers equation: E - e*sinE = Me
%
% Inputs: 
%   Me    Mean anomaly (rad)
%   e     Eccentricity (0 <= e < 1)
%   tol   Error tolerance. Optional.
%
% Outputs:
%   E     Eccentric anomaly (rad)
%

if nargin<3
  tol = 1e-8;
end

% initial guess
if( Me>pi )
  E = Me - e/2;
else
  E = Me + e/2;
end

dist = 1;
while dist>tol
  
  % compute function value and function derivative value
  f = E-e.*sin(E)-Me;
  fp = 1-e.*cos(E);

  % distance between current and next iterate
  dist = max(abs(f./fp));
  
  E = E-f./fp;

end


