function [x,y] = ellipse_locus(F1,F2,a)

% find the locus of points for an ellipse with SMA a and two foci locations
% at F1, F2

if nargin<1 
  a = 1;
  e = 0.1;
  F1 = [0;0];
  F2 = [-2*a*e;0];
end

% constraint is r1 + r2 = 2*a
% where r1 and r2 are the distances from F1 and F2 to the point,
% respectively

