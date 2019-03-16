function [f,g,fdot,gdot] = LagrangeCoeffZ( r1m, r2m, dTheta, z, mu )

% Compute the Lagrange coefficients in terms of universal variable z
%
%   Inputs: 
%     r1m     Magnitude of position vector r1  
%     r2m     Magnitude of position vector r2
%     dTheta  Angle between vector r1 and r2
%     z       Universal variable "z"
%     mu      Gravitational constant 
% 
% 	Outputs:
%     f   
%     g
%     fdot
%     gdot 
%


C     = stumpC(z);
S     = stumpS(z);
A     = sin(dTheta)*sqrt(r1m*r2m/(1-cos(dTheta)));  % 5.35
y     = r1m+r2m+A*( z*S-1 )/sqrt(C);                % 5.38
f     = 1-y/r1m;                                    % 5.46a
g     = A*sqrt(y/mu);                               % 5.46b
fdot  = sqrt(mu)/(r1m*r2m)*sqrt(y/C)*(z*S-1);       % 5.46c
gdot  = 1-y/r2m;                                     % 5.46d