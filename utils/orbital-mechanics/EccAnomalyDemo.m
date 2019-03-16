%% How does E vary with True anomaly ... 
% E: Eccentric anomaly for 0 < e < 1

% Define a range of true anomaly "theta"
theta = linspace(0,2*pi);

% define a few eccentricity values
eSet = [0.1 0.3 0.5 0.7 0.9];

% size of our vectors
nTh = length(theta);
nEcc = length(eSet);

% initialize a matrix with "nEcc" rows and "nTh" columns
E = zeros(nEcc,nTh);

% compute the eccentric anomaly
for i=1:nEcc
  E(i,:) = EccAnomFromTrueAnom(theta, eSet(i) );
end
E( E<0 ) = E( E<0 )+2*pi;

%% Plot eccentric anomaly
figure('name','Eccentric Anomaly vs. True Anomaly')
plot(theta,E)
legend(strread(num2str(eSet),'%s'))
grid on
set(gca,'fontsize',14)
xlabel('True Anom (rad)')
ylabel('Ecc Anom (rad)')

%% How does E vary with Mean anomaly ... 
% E: Eccentric anomaly for 0 < e < 1

% Define a range of mean "Me"
Me = linspace(0,2*pi);

% define a few eccentricity values
eSet = [0.1 0.3 0.5 0.7 0.9];

% size of our vectors
nM = length(Me);
nEcc = length(eSet);

% initialize a matrix with "nEcc" rows and "nM" columns
E = zeros(nEcc,nM);

% compute the eccentric anomaly
for i=1:nEcc
  for j=1:nM
      E(i,j) = EccAnomFromMeanAnom(Me(j), eSet(i) );
  end
end
E( E<0 ) = E( E<0 )+2*pi;

% Plot eccentric anomaly
figure('name','Eccentric Anomaly vs. Me')
plot(Me,E)
legend(strread(num2str(eSet),'%s'))
grid on
set(gca,'fontsize',14)
xlabel('Mean Anom (rad)')
ylabel('Ecc Anom (rad)')


%% How does F vary with True anomaly ... 
% F: HYPERBOLIC Eccentric Anomaly for e > 1

% pick an e > 1
e = 1.5;

% define a true anomaly range
theta = linspace(2.2, 2.4, 1e7); % is this a good idea? why not?

% compute 
F = HypEccAnomFromTrueAnom( theta, e );

figure('name','Hyperbolic Eccentric Anomaly vs. True Anomaly')
plot(theta,F)
grid on
set(gca,'fontsize',14)
xlabel('True Anom. (rad)')
ylabel('Hyp. Ecc. Anom. (rad)')

% we need to enforce: theta < theta_inf
theta_inf = acos( -1/e );

% draw it on the plot
hold on
plot([theta_inf theta_inf],[0,max(F)],'k--','linewidth',3)

%% How does F vary with Mean anomaly ... 
% F: HYPERBOLIC Eccentric Anomaly for e > 1

% pick an e > 1
e = 3.5;

% Let's plot Mh vs F
F = linspace(0,2*pi);
Mh = e*sinh(F) - F;

%figure('name','Hyp. Mean Anomaly vs. Hyp. Ecc. Anomaly')
semilogy(F,Mh,'linewidth',3) % plot with Mh on log scale
grid on, hold on
set(gca,'fontsize',14)
ylabel('Mean Anom. (rad)')
xlabel('Hyp. Ecc. Anom. (rad)')


