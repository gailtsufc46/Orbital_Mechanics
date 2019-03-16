function Q = GeoEqToPerifocal( inc, RA, w )

% Compute the matrix that rotates from the geocentric equatorial frame to
% the perifocal frame
%
% Inputs:
%   inc   Inlination angle (rad)
%   RA    Right ascension of the ascending node (rad)
%   w     Argument of perigee (rad)
%

Q = EulerRot(3,w)*EulerRot(1,inc)*EulerRot(3,RA);
