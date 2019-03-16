%% Illustrate how F changes with theta and M

% define an eccentricy (e>1) and a semi-major axis
mu = 398600.4;
e = 2;
a = 8000;

% compute a vector of true anomaly values between +/- theta_inf
thinf = acos(-1/e);
thx = linspace(-thinf+0.01, thinf-0.01, 1e6);

% compute the hyperbolic eccentric anomaly values
Fx = HypEccAnomFromTrueAnom( thx, e );

% compute the hyperbolic mean anomaly values
Mx = e*sinh(Fx)-Fx;

% compute time (need mean orbit rate, "n")
n = sqrt(mu/a);
tx = Mx/n;

% plot F vs. Mh
figure
plot(Mx*180/pi,Fx*180/pi,'linewidth',2)
set(gca,'fontsize',14)
xlabel('Hyp Mean Anom (deg)')
ylabel('Hyp Ecc Anom (deg)')
grid on

% plot F and Mh together, vs. theta
figure
plot(thx*180/pi,Fx*180/pi,...
  thx*180/pi,Mx*180/pi,'linewidth',2)
set(gca,'fontsize',14)
xlabel('True Anom (deg)')
ylabel('F(\theta) and M_h(\theta)')
grid on
legend('Hyp Ecc Anom, F','Hyp Mean Anom, M_h')
axis([-100 130 -300 1200])
set(gca,'xtick',[-90 : 30 : 120])

% plot F, Mh and theta vs. time
figure
plot(tx,Mx,tx,Fx,tx,thx), legend('Mean anomaly (M)','Hyp. Ecc. Anomaly (F)','True Anomaly')
