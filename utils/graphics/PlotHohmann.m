function d = PlotHohmann( rA, rF, mu )

if( nargin<1 )
  rA = 7000;
  rF = 12000;
  mu = 398600.44;
  d = PlotHohmann( rA, rF, mu );
  return
end

q = linspace(0,2*pi);

figure
plot(rA*cos(q),rA*sin(q),'b','linewidth',2)
hold on
plot(rF*cos(q),rF*sin(q),'r','linewidth',2)
axis equal
grid on
xlabel('x')
ylabel('y')
plot(0,0,'k.','markersize',30)

g = OrbMvrHohmann(rA,rF,rA,rF,mu);

[rT,xT,yT] = orbital_path( g.a, g.e, linspace(0,pi) );
plot(xT,yT,'k--','linewidth',2)

