function d = PlotBiElliptic( rA, rB, rC, mu )

if( nargin<1 )
  rA = 7000;
  rB = 30000;
  rC = 12000;
  mu = 398600.44;
  d = PlotBiElliptic( rA, rB, rC, mu );
  return
end

q = linspace(0,2*pi);

figure
plot(rA*cos(q),rA*sin(q),'b','linewidth',2)
hold on
plot(rC*cos(q),rC*sin(q),'r','linewidth',2)
axis equal
grid on
xlabel('x')
ylabel('y')
plot(0,0,'k.','markersize',30)

g = OrbMvrBiElliptic(rA,rB,rC,mu);

[rT,xT,yT] = orbital_path( g.a2, g.e2, linspace(0,pi) );
plot(xT,yT,'k--','linewidth',2)

[rT,xT,yT] = orbital_path( g.a3, g.e3, linspace(pi,2*pi) );
plot(xT,yT,'m--','linewidth',2)
