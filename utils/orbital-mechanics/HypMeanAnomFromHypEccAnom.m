function Mh = HypMeanAnomFromHypEccAnom( F, e )
%
% Compute hyperbolic mean anomaly for a given hyperbolic eccentric anomaly
%
% Inputs: 
%   F         Hyperbolic eccentric anomaly
%   e         Eccentricity (e > 1)
%
% Outputs:
%   Mh        Hyperbolic ean anomaly
%

Mh = e*sinh(F) - F;
