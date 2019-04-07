function [v1, v2, dTheta] = solveLambert(r1,r2,tof,mu,dir)
% Written by Garrett Ailts
%
% Usage: [v1, v2, dTheta] = solveLambert(r1,r2,tof,mu,dir)
%
% Description: Function that solves Lambert's Problem given two points in 
% space and a time of flight. Returns the orbital elements for a Keplerian 
% orbit that crosses r1 and r2 in the given time of flight tof. The
% returned orbit can be for the prograde or retrograde solution and can be
% specified with the dir input. 
%
% Inputs: r1 - position vector specifying the first point in space (km)
%         r2 - position vector specifying the second point in space (km)
%        tof - time of flight between r1 and r2 (s)
%         mu - gravitational parameter of central body (km^3/s^2)
%        dir - string specifying the direction of the orbit ('prograde' or
%        'retrograde')
%
% Output: a - semi major axis (km)
%         e - eccentricity
%     theta - true anomaly (rad)
%     OMEGA - right ascension of ascending node (rad)
%     omega - arguement of perigee (rad)
%       inc - inclination (rad)
%         h - specific angualr momentum (km^2/s)

if nargin < 5
    dir = 'prograde';
end
if nargin < 4
    mu = 398600.44;
end
if tof <= 0
    error("Time of flight must be greater than 0!");
end

%% Compute Vector Magnitudes and Delta Theta
r1mag = norm(r1);
r2mag = norm(r2);
n = cross(r1,r2);
dTheta = acos(dot(r1,r2)/r1mag/r2mag);
if strcmp(dir,'prograde')==1
    if (n(3)<0)
       dTheta = 2*pi-dTheta;
    end
elseif strcmp(dir,'retrograde')
    if(n(3)>=0)
        dTheta = 2*pi-dTheta;
    end
end
        

%% Create Function Handles For Constituents of Universal Kepler's Eq
A = sin(dTheta)*sqrt(r1mag*r2mag/(1-cos(dTheta)));
C = @(z) stumpffC(z);
S = @(z) stumpffS(z);
y = @(z) r1mag+r2mag+A*(z.*S(z)-1)./sqrt(C(z));

%% Define Universal Kepler's Eq and it's Derivative.
F = @(z) (y(z)./C(z)).^(3/2).*S(z)+A*sqrt(y(z))-sqrt(mu)*tof;
Fprime = @(z) DFofZ(r1mag, r2mag, A, z);
           
%% Find Initial Guess for z and Calculate z
z0s = FindInitialZGuessForLambert(y,F);
for i = 1:length(z0s)
    [z,~,success] = NewtRaph(F,Fprime,z0s(i),1e-8,1000);
    if success == 1
        break;
    end
end

%% Compute Neccessary Lagrange Coeffiecients
[f, g, ~, gdot] = lagrangeFromZ(r1mag, r2mag, A, z, mu);

%% Compute v1 and v2
v1 = (1/g)*(r2-f*r1);
v2 = (1/g)*(gdot*r2-r1);
end


        
        