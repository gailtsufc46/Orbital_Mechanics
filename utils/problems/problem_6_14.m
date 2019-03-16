%% Problem 6.14 from the textbook
% Sequence of Hohmann transfers in LEO for the space shuttle

Re = 6378.14;
mu = 398600.44;

d12a = OrbMvrHohmann(Re+302,Re+291,Re+296,Re+259,mu);
d12b = OrbMvrHohmann(Re+302,Re+259,Re+296,Re+291,mu);

d23 = OrbMvrHohmann(Re+291,Re+259,Re+259,Re+259,mu);

d34 = OrbMvrHohmann(Re+259,Re+255,Re+259,Re+194,mu);

dVTot = ...
  min([d12a.dVTot,d12a.dVTotx,d12b.dVTot,d12b.dVTotx]) + ...
  min([d23.dVTot, d23.dVTotx]) + ...
  min([d34.dVTot, d34.dVTotx]);

minTotDV = dVTot*1e3
