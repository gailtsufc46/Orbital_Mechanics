function d = HohmannDV( r1, r2, mu )

%   HohmannDV - Compute the delta-vs for a Hohmann transfer
% 
%   Inputs:
%     r1      Orbit 1 radius (km)
%     r2      Orbit 2 radius (km)
%     mu      Gravitational constant of the central body (km^3/s^2)
%
%   Outputs:
%     d       Data structure with fields:
%                v1: Orbit 1 velocity
%                v2: Orbit 2 velocity
%                hT: Transfer orbit angular momentum (specific)
%               vT1: Transfer orbit velocity at position 1
%               vT2: Transfer orbit velocity at position 2
%               dV1: Delta-v at position 1
%               dV2: Delta-v at position 2
%



if( nargin < 1 )
    r1 = 6378+480;
    r2 = 6478+16000;
    mu = 398600.44;
end

d.v1 = sqrt(mu/r1);
d.v2 = sqrt(mu/r2);

d.hT = sqrt(2*mu)*sqrt( r1*r2/(r1+r2) );

d.vT1 = sqrt(2*mu)*sqrt( (r2/r1) / (r2+r1) );
d.vT2 = sqrt(2*mu)*sqrt( (r1/r2) / (r2+r1) );

d.dV1 = d.vT1-d.v1;
d.dV2 = d.v2-d.vT2;
