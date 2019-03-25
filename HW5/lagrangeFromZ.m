function [f, g, fdot, gdot] = lagrangeFromZ(r1mag, r2mag, A, z, mu)
% Written by Garrett Ailts
%
% Usage: [f, g, fdot, gdot] = lagrangeFromZ(r1mag, r2mag, A, yz);
% 
% Inputs: r1mag - magnitude of r1
%         r2mag - magnitude of r2
%             A - constant of universal variable equation
%            yz - y(z) evaluated at z
%
% Outputs: f - Lagrange Coefficient f
%          g - Lagrange Coefficient g
%       fdot - time derivative of f
%       gdot - time derivative of g

%% Evaluate Stumpff Functions and y(z) @z
C = stumpffC(z);
S = stumpffS(z);
yz = r1mag+r2mag+A*(z*S-1)/sqrt(C);


%% Compute Lagrange Coefficients
f = 1-yz/r1mag;
g = A*sqrt(yz/mu);
fdot = sqrt(mu)*sqrt(yz/C)*(z*S-1)/r1mag/r2mag;
gdot = 1-yz/r2mag;