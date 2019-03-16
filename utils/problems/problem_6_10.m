%% Problem 6.10 from the textbook
% Hohmann transfer from Earth to Mars

muSun = 1.327e11; % km^3/sec^2
rE    = 1.496e8;  % km
rM    = 2.279e8;  % km

d = OrbMvrHohmann( rE, rM, rE, rM, muSun );

time = d.T/2;
timeDays = time/86400

nMars = OrbRate(rM,muSun);

alpha = pi - nMars*time;
alphaDeg = alpha*180/pi