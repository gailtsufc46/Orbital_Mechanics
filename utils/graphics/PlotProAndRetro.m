function PlotProAndRetro( r1, r2, v1, v1r, mu )

% Plot the PRO- and RETRO-grade orbits together

[xI,yI,zI] = GenerateOrbitPositionFromRV( r1, v1, mu );
[xIR,yIR,zIR] = GenerateOrbitPositionFromRV( r1, v1r, mu );

figure, 
plot3(0,0,0,'k.','markersize',20,'displayname','Central Body'), hold on, grid on
plot3(xI,yI,zI,'b','linewidth',2,'displayname','Prograde Orbit'), 
plot3(xIR,yIR,zIR,'r','linewidth',2,'displayname','Retrograde Orbit'), 
plot3(r1(1),r1(2),r1(3),'g.','markersize',20,'displayname','R1')
plot3(r2(1),r2(2),r2(3),'m.','markersize',20,'displayname','R2')
axis equal
xlabel('x_{ECI} (km)')
ylabel('y_{ECI} (km)')
zlabel('z_{ECI} (km)')
legend
