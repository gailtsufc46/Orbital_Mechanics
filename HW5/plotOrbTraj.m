function plotOrbTraj(COE,thetaVec,mu)
% Written by Garrett Ailts
%
% Usage: plotOrbTraj(COE,thetaVec)
%
% Description: Function generates a 3D plot of an orbit given common
% orbital elements and a vector of true anomalies.
%
% Inputs: COE - Structure containing common orbital elements
%            Contents of COE:
%                       a - semi major axis (km)
%                       e - eccentricity
%                   theta - true anomaly (rad)
%                   OMEGA - right ascension of ascending node (rad)
%                   omega - arguement of perigee (rad)
%                     inc - inclination (rad)
%
% Outputs: 3D figure showing the orbital trajectory (km)

%% Preallocation
r = zeros(length(thetaVec),3);
%% Convert COE and True Anomaly Vector Into Inertial Coordinates
for i = 1:length(thetaVec)
    [r(i,:), ~, ~] = orbEl2rv(COE.a, COE.e, thetaVec(i), COE.OMEGA, COE.omega, COE.inc, mu);
end

%% Plot Orbital Trajectory
plot3(r(:,1),r(:,2),r(:,3))
