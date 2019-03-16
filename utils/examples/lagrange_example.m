
factor = 1.5;

r0 = [7000; 100; 0]; 
v0 = [4; 8; 0] * factor;
mu = 398600.44;   % grav. param. for Earth (km^3/s^2)

dTh = linspace(0,2*pi);

r = zeros(3,length(dTh));
v = r;

for i=1:length(dTh)
  
  [f,g,fdot,gdot] = LagrangeCoeff(r0,v0,mu,dTh(i));
  
  r(:,i) = f*r0 + g*v0;
  v(:,i) = fdot*r0 + gdot*v0;
  
end

figure, plot(r(1,:),r(2,:)), axis equal, grid on
