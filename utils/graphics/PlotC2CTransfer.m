function [vT1,vT2,dv1,dv2] = PlotC2CTransfer( rA, rB, thA, thB, t0, t1, t2, mu, name )

% Plot a transfer from orbit A to orbit B.
% The orbits have radii rA and rB.
% Orbiting bodies are initially at the given phase angles, thA and thB.
% Wait until time t1.
% The transfer starts at time t1 and ends at time t2.
% mu is the gravitational parameter of the central body.

if( nargin<9 )
  name = 'Circle to Circle Transfer';
end

% prevent numerical issues when thA-thB = PI
if( abs(thA-thB)-pi < 1e-8 )
  thB = thB+1e-4;
end

thx = linspace(0,2*pi,300);

figure('name',name)

plot(0,0,'y.','markersize',30,'handlevisibility','off'), hold on, grid on, zoom on
axis equal

plot(rA*cos(thx),rA*sin(thx),'b','handlevisibility','off')
plot(rB*cos(thx),rB*sin(thx),'r','handlevisibility','off')

plot(rA*cos(thA),rA*sin(thA),'bo','markersize',14,'displayname','A at t0')
plot(rB*cos(thB),rB*sin(thB),'ro','markersize',14,'displayname','B at t0')

nA = sqrt(mu/rA^3);
nB = sqrt(mu/rB^3);

thA1 = thA+nA*(t1-t0);
thA2 = thA+nA*(t2-t0);
thB1 = thB+nB*(t1-t0);
thB2 = thB+nB*(t2-t0);

plot(rA*cos(thA1),rA*sin(thA1),'bo','markersize',14,'handlevisibility','off')
plot(rB*cos(thB1),rB*sin(thB1),'ro','markersize',14,'handlevisibility','off')
plot(rA*cos(thA1),rA*sin(thA1),'bx','markersize',14,'displayname','A at t1')
plot(rB*cos(thB1),rB*sin(thB1),'rx','markersize',14,'displayname','B at t1')
plot(rA*cos(thA2),rA*sin(thA2),'b.','markersize',32,'displayname','B at t1')
plot(rB*cos(thB2),rB*sin(thB2),'r.','markersize',32,'displayname','B at t2')

rA1 = rA*[cos(thA1);sin(thA1);0];
rB2 = rB*[cos(thB2);sin(thB2);0];

vA1 = sqrt(mu/rA)*[-sin(thA1);cos(thA1);0];
vB2 = sqrt(mu/rB)*[-sin(thB2);cos(thB2);0];

[vT1,vT2]=LambertSolver(rA1,rB2,t2-t1,mu);
dv1 = vT1 - vA1;
dv2 = vT2 - vB2;

tt = linspace(0,t2-t1,250);
rTrans = RVAtTFromR0V0(rA1,vT1,tt,mu);
plot(rTrans(1,:),rTrans(2,:),'m','displayname','Transfer orbit')
%legend

text(rA*cos(thA)*1.1,rA*sin(thA)*1.1,sprintf('%2.1f days',t0/86400))
text(rB*cos(thB)*1.1,rB*sin(thB)*1.1,sprintf('%2.1f days',t0/86400))
text(rA*cos(thA1)*1.1,rA*sin(thA1)*1.1,sprintf('%2.1f days',t1/86400))
text(rB*cos(thB1)*1.1,rB*sin(thB1)*1.1,sprintf('%2.1f days',t1/86400))
text(rA*cos(thA2)*1.1,rA*sin(thA2)*1.1,sprintf('%2.1f days',t2/86400))
text(rB*cos(thB2)*1.1,rB*sin(thB2)*1.1,sprintf('%2.1f days',t2/86400))
