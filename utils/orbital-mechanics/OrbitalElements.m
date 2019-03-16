function [a,inc,W,w,e,th] = OrbitalElements( el )

% Extract the individual elements from a vector of orbital elements.
%
% Inputs:
%   el            Orbital elements vector [a,i,W,w,e,th]
%
% Outputs:
%   a             Semi major axis     [km]
%   inc           Inclination         [rad]
%   W             Right ascension     [rad]
%   w             Argument of perigee [rad]
%   e             Eccentricity 
%   th            True anomaly        [rad]
%


a = el(1);
inc = el(2);
W = el(3);
w = el(4);
e = el(5);
th = el(6);
