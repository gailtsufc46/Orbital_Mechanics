function [dr,dv] = CWEqns( dr0, dv0, n, t )

%   Clohessy-Wiltshire equations. Compute the relative position and velocity
%   over time given the initial relative position and velocity, and the mean
%   orbit rate of the reference orbit. This method is valid only for circular
%   reference orbits, and small deviations from that orbit.
% 
%   The relative motion is defined in a rotating coordinate frame attached to
%   the reference satellite. This is referred to as a local-vertical,
%   local-horizontal (LVLH) frame. 
% 
%   The x,y,z directions in LVLH are defined as follows:
%     x: Radial
%     y: Along-Track
%     z: Cross-Track
% 
%   Inputs:
%     dr0     Initial relative position vector in LVLH  
%     dr0     Initial relative velocity vector in LVLH  
%     n       Mean orbit rate of reference orbit        (rad/sec)
%     t       Time elapsed since initial state          (sec)

x0 = dr0(1);
y0 = dr0(2);
z0 = dr0(3);
u0 = dv0(1);
v0 = dv0(2);
w0 = dv0(3);

nt = n*t;

c = cos(nt);
s = sin(nt);

dx = (4-3*c)*x0 + s/n*u0 + 2/n*(1-c)*v0;
dy = 6*(s-nt)*x0 + y0 + 2/n*(c-1)*u0 + 1/n*(4*s-3*nt)*v0;
dz = c*z0 + 1/n*s*w0;

du = 3*n*s*x0 + c*u0 + 2*s*v0;
dv = 6*n*(c-1)*x0 - 2*s*u0 + (4*c-3)*v0;
dw = -n*s*z0 + w0*c;

dr = [dx;dy;dz];
dv = [du;dv;dw];
