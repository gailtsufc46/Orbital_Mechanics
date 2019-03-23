% HW5 Script
% Orbital Mechanics
% Garrett Ailts

%% Constants
mu = 398600.44;
%% Define Cases
% Types 
dir = ['prograde','retrograde'];

% Case 1
cases(1).r1 = [-3472.57238396015;8042.69586228059;2196.56246176354];
cases(1).r2 = [-8710.68403512538;-3815.7502616566;885.874473990578];
cases(1).dT = 2157.063791;
% Case 2
cases(2).r1 = [-5192.64707332683;-2730.67520032545;-1312.84187347459];
cases(2).r2 = [-1396.17652892469;-2398.88103277394;7822.986845429];
cases(2).dT = 1288.841280;
% Case 3
cases(3).r1 = [506.717068361234;-470.869269481925;138.096930726997];
cases(3).r2 = [-10298.3068750911;6203.57304746482;-16637.5263693708];
cases(3).dT = 1800.000000;

%% Solve Lamberts for Each Case
for i = 1:length(cases)
    for j = 1:length(dir)
    [v1, v2, dTheta] = solveLambert(cases(i).r1,cases(i).r2,cases(i).dT,mu,dir(j));
    cases(i).dir(j).v1 = v1; % dir(1) is prograde, dir(2) is retrograde
    cases(i).dir(j).v2 = v2;
    cases(i).dir(j).dTheta =dTheta;
    end
end

%% Obtain COE and Initial and Final True Anomaly for Each Case
for i = 1:length(cases)
    for j = 1:length(dir)
    [a, e, theta, OMEGA, omega, inc, h] = rv2OrbEl(cases(i).r1,cases(i).dir(j).v1,mu);
    cases(i).dir(j).a = a;
    cases(i).dir(j).e =e ;
    cases(i).dir(j).theta = theta;
    cases(i).dir(j).OMEGA = OMEGA;
    cases(i).dir(j).omega = omega;
    cases(i).dir(j).inc = inc;
    cases(i).dir(j).h = h;
    cases(i).dir(j).thetaf = theta+cases(i).dir(j).dTheta;
    end
end

%% Plot Lambert Solutions
thetaVec = linspace(0,2*pi,100);
for i = 1:length(cases)
    figure(i), hold on
    for j = 1:length(dir)
        plotOrbTraj(cases(i).dir(j), thetaVec, mu);
    end
    xlabel('I (km)'), ylabel('J (km)'), zlabel('K (km)');
    titel(sprintf("Lambert's Solutions For Case %u",i));
    legend('Prograde', 'Retrograde'), hold off
end
    








