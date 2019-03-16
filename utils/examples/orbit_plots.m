%% Generate some sample orbit plots

[r,x,y]    = orbital_path(1,0,linspace(0,2*pi));
[r2,x2,y2] = orbital_path(1,0.5,linspace(0,2*pi));
[r3,x3,y3] = orbital_path(1,0.75,linspace(0,2*pi));
[r4,x4,y4] = orbital_path(1,0.9,linspace(0,2*pi,1e3));
figure, 
plot(x,y), axis equal, zoom on, grid on, hold on
plot(x2,y2,'r')
plot(x3,y3,'c')
plot(x4,y4,'m')
plot(0,0,'k.','markersize',28)

th1 = linspace(0,2*pi,1e3);
th2 = linspace(-.8*pi,.8*pi);
th3 = linspace(-.6*pi,.6*pi);
a = 1;
e = 0.0; p = a*(1-e^2); [rc,xc,yc] = orbital_path(p,e,th1);
e = 0.5; p = a*(1-e^2); [re,xe,ye] = orbital_path(p,e,th1);
e = 1.0; p = 0.5;       [rp,xp,yp] = orbital_path(p,e,th2);
e = 2.0; p = 0.5;       [rh,xh,yh] = orbital_path(p,e,th3);
figure
plot(xc,yc,xe,ye,xp,yp,xh,yh,'linewidth',2)
axis equal, zoom on, grid on, hold on
plot(0,0,'k.','markersize',28)
set(gca,'fontsize',14)
legend('e = 0, Circle','e = 0.5, Ellipse','e = 1, Parabola','e = 2, Hyperbola')
