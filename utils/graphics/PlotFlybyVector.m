function PlotFlybyVector( f )

% Plot the velocity vectors associated with a flyby maneuver
% 
% Inputs:
%   f   Data structure output from Flyby.m
%
% Outputs:
%   none
%

figure
% plot the planet heliocentric velocity
plot([0 f.vPH(1)],[0 f.vPH(2)],'k','linewidth',2)
hold on
grid on
axis equal

% plot a circle around the tip of the planet velocity, radius = "vinf"
ang = linspace(0,2*pi);
xc = f.vinf*cos(ang);
yc = f.vinf*sin(ang);
plot(xc+f.vPH(1),yc+f.vPH(2),'r')

% plot the initial spacecraft velocity
plot([0 f.v0H(1)],[0 f.v0H(2)],'b--')

% plot the final spacecraft velocity
plot([0 f.vFH(1)],[0 f.vFH(2)],'b')

xlabel('x (helio)')
ylabel('y (helio)')

legend('Planet velocity','Locus of possible spacecraft velocities',...
  'Initial spacecraft velocity','Final spacecraft velocity after turn',...
  'location','best')
